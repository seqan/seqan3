// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::coordinate_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/view_iota_simd.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

namespace seqan3::detail
{

//------------------------------------------------------------------------------
// coordinate_matrix
//------------------------------------------------------------------------------

/*!\brief A matrix over coordinates.
 * \ingroup alignment_matrix
 * \implements std::ranges::forward_range
 *
 * \tparam index_t The underlying matrix index type; must mode std::integral or seqan3::simd::simd_concept.
 *
 * \details
 *
 * This matrix emulates a two-dimensional index, such that each cell inside of the alignment matrix can be located
 * through a unique coordinate (=index). In the alignment algorithm this matrix is paired with the alignment matrix
 * to form an indexed alignment matrix.
 *
 * This matrix is cheap as it stores only the dimensions of the matrix and does not allocate any memory for the
 * coordinates. It uses the seqan3::detail::matrix_coordinate_iterator to iterate
 * over the columns of this virtual matrix. When the seqan3::detail::coordinate_matrix::matrix_coordinate_iterator
 * is dereferenced it returns an on-the-fly constructed range
 * ([prvalue](https://en.cppreference.com/w/cpp/language/value_category))
 * representing the current indexed column. The reference type of this range is seqan3::detail::matrix_coordinate.
 *
 * ### Simd mode
 *
 * If the `index_t` is a simd vector the matrix implements a vectorised coordinate matrix instead. In this case the
 * matrix uses a seqan3::detail::coordinate_matrix::matrix_coordinate_simd_iterator to iterate over the vectorised
 * matrix. When the seqan3::detail::coordinate_matrix::matrix_coordinate_iterator is dereferenced it returns an
 * on-the-fly constructed range ([prvalue](https://en.cppreference.com/w/cpp/language/value_category)) representing the
 * current simd coordinate column. The reference type of this range is seqan3::detail::simd_matrix_coordinate.
 */
template <typename index_t>
//!\cond
    requires std::integral<index_t> || simd_concept<index_t>
//!\endcond
class coordinate_matrix
{
private:

    /*!\brief A function object that converts a column index and a row index to a seqan3::detail::matrix_coordinate.
     * \details
     *
     * This helper function object is used for the seqan3::detail::coordinate_matrix::coordinate_matrix_iterator and
     * seqan3::detail::coordinate_matrix::coordinate_matrix_simd_iterator to convert the column and row index into a
     * seqan3::detail::matrix_index.
     */
    struct convert_to_matrix_coordinate
    {
        //!\brief The index of the represented column.
        index_t column_index{};

        /*!\brief The conversion operator
         *
         * \param[in] row_index The index of the represented row.
         *
         * \returns seqan3::detail::matrix_index for the current column and row index.
         */
        auto operator()(index_t const row_index) noexcept
        {
            return matrix_index{row_index_type{row_index}, column_index_type{column_index}};
        }
    };

    // The coordinate matrix iterator for simd index types.
    template <simd_concept simd_index_t>
    class coordinate_matrix_simd_iterator;
    // The coordinate matrix iterator for regular integral types.
    class coordinate_matrix_iterator;

    //!\brief Helper type alias to conditionally extract the scalar type from the seqan3::simd::simd_traits type.
    template <typename lazy_simd_traits_t>
    using lazy_scalar_type_t = typename instantiate_t<lazy_simd_traits_t>::scalar_type;

    //!\brief The internal size type which depends on `index_t` being a simd vector or a scalar type.
    using size_type = lazy_conditional_t<simd_concept<index_t>,
                                         lazy<lazy_scalar_type_t, lazy<simd_traits, index_t>>,
                                         size_t>;

    //!\brief The iterator type which depends on `index_t` being a simd vector or a scalar type.
    using iterator_type = lazy_conditional_t<simd_concept<index_t>,
                                             lazy<coordinate_matrix_simd_iterator, index_t>,
                                             coordinate_matrix_iterator>;

    //!\brief The number of columns (corresponds to the size of one row).
    size_type column_count{};
    //!\brief The number od rows (corresponds to the size of one column).
    size_type row_count{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    coordinate_matrix() = default; //!< Defaulted.
    coordinate_matrix(coordinate_matrix const &) = default; //!< Defaulted.
    coordinate_matrix(coordinate_matrix &&) = default; //!< Defaulted.
    coordinate_matrix & operator=(coordinate_matrix const &) = default; //!< Defaulted.
    coordinate_matrix & operator=(coordinate_matrix &&) = default; //!< Defaulted.
    ~coordinate_matrix() = default; //!< Defaulted.

    /*!\brief Resets the coordinate matrix with the given end column index and end row index representing the new
     *        dimensions of the matrix.
     *\tparam row_count_t The size type for the column count; must model std:integral.
     *\tparam row_size_t The size type for the row count; must model std:integral.
     *
     * \param[in] column_count \copybrief seqan3::detail::coordinate_matrix::column_count
     * \param[in] row_count \copybrief seqan3::detail::coordinate_matrix::row_count
     *
     * \details
     *
     * ### Complexity
     *
     * Constant
     *
     * ### Exception
     *
     * noexcept
     */
    template <std::integral row_count_t, std::integral row_size_t>
    void resize(column_index_type<row_count_t> const column_count,
                row_index_type<row_size_t> const row_count) noexcept
    {
        this->column_count = column_count.get();
        this->row_count = row_count.get();
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column of the matrix.
    iterator_type begin() const noexcept
    {
        return iterator_type{static_cast<size_type>(0), row_count};
    }

    //!\brief Returns the iterator pointing to the end column of the matrix.
    iterator_type end() const noexcept
    {
        return iterator_type{column_count, row_count};
    }
    //!\}
};

//------------------------------------------------------------------------------
// coordinate_matrix_iterator
//------------------------------------------------------------------------------

/*!\brief The iterator for the seqan3::detail::coordinate_matrix.
 * \implements std::forward_iterator
 *
 * \details
 *
 * Iterates over the columns of the underlying seqan3::detail::coordinate_matrix. The iterator returns a transformed
 * std::ranges::views::iota view over the row indices. In the transformation every row index is converted to a
 * seqan3::detail::matrix_coordinate using the seqan3::detail::convert_to_matrix_coordinate function object.
 */
template <typename index_t>
//!\cond
    requires std::integral<index_t> || simd_concept<index_t>
//!\endcond
class coordinate_matrix<index_t>::coordinate_matrix_iterator
{
private:
    //!\brief The currently represented column index.
    size_t column_id{0};
    //!\brief The number of rows (corresponds to the size of a column).
    size_t row_count{0};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::views::iota(size_t{0u}, row_count)
                              | std::views::transform(convert_to_matrix_coordinate{column_id}));
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::forward_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    coordinate_matrix_iterator() noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator(coordinate_matrix_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator(coordinate_matrix_iterator &&) noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator & operator=(coordinate_matrix_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_iterator & operator=(coordinate_matrix_iterator &&) noexcept = default; //!< Defaulted.
    ~coordinate_matrix_iterator() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the current column index and the row index marking the end
     *        of the rows (size of one column).
     *
     * \param[in] column_id \copybrief seqan3::detail::coordinate_matrix_iterator::column_id
     * \param[in] row_count \copybrief seqan3::detail::coordinate_matrix_iterator::row_count
     */
    explicit coordinate_matrix_iterator(size_t const column_id, size_t const row_count) noexcept :
        column_id{column_id},
        row_count{row_count}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    auto operator*() const
    {
        return std::views::iota(size_t{0u}, row_count)
             | std::views::transform(convert_to_matrix_coordinate{column_id});
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator to the next column.
    coordinate_matrix_iterator & operator++()
    {
        ++column_id;
        return *this;
    }

    //!\brief Increments the iterator to the next column and returns the iterator pointing to the previous column.
    coordinate_matrix_iterator operator++(int)
    {
        coordinate_matrix_iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(coordinate_matrix_iterator const & lhs, coordinate_matrix_iterator const & rhs)
    {
        return lhs.column_id == rhs.column_id;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(coordinate_matrix_iterator const & lhs, coordinate_matrix_iterator const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}
};

//------------------------------------------------------------------------------
// coordinate_matrix_simd_iterator
//------------------------------------------------------------------------------

/*!\brief The iterator for the vectorised seqan3::detail::coordinate_matrix.
 * \implements std::forward_iterator
 *
 * \tparam index_t The type of the index; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * Iterates over the columns of the underlying seqan3::detail::coordinate_matrix. The iterator returns a transformed
 * seqan3::views::iota_simd view over the row indices. In the transformation every row index is converted to a
 * seqan3::detail::matrix_index using the seqan3::detail::convert_to_matrix_coordinate function object.
 */
template <typename index_t>
//!\cond
    requires std::integral<index_t> || simd_concept<index_t>
//!\endcond
template <simd_concept simd_index_t>
class coordinate_matrix<index_t>::coordinate_matrix_simd_iterator
{
private:
    //!\brief The scalar type.
    using scalar_type = typename simd_traits<simd_index_t>::scalar_type;

    //!\brief The column index as simd vector.
    simd_index_t simd_column_index{};
    //!\brief The currently represented column index.
    scalar_type scalar_column_index{};
    //!\brief The end index of the row.
    scalar_type column_size{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(views::iota_simd<simd_index_t>{static_cast<scalar_type>(0), column_size}
                              | std::views::transform(convert_to_matrix_coordinate{}));
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::forward_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    coordinate_matrix_simd_iterator() noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator(coordinate_matrix_simd_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator(coordinate_matrix_simd_iterator &&) noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator & operator=(coordinate_matrix_simd_iterator const &) noexcept = default; //!< Defaulted.
    coordinate_matrix_simd_iterator & operator=(coordinate_matrix_simd_iterator &&) noexcept = default; //!< Defaulted.
    ~coordinate_matrix_simd_iterator() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the current column index and the row index marking the end
     *        of the rows (size of one column).
     *
     * \param[in] column_index The column index to point to.
     * \param[in] column_size The size of the column.
     */
    explicit coordinate_matrix_simd_iterator(scalar_type const column_index, scalar_type const column_size) noexcept :
        simd_column_index{simd::fill<simd_index_t>(column_index)},
        scalar_column_index{column_index},
        column_size{column_size}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    reference operator*() const
    {
        return views::iota_simd<simd_index_t>{static_cast<scalar_type>(0), column_size}
             | std::views::transform(convert_to_matrix_coordinate{simd_column_index});
    }
    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator to the next column.
    coordinate_matrix_simd_iterator & operator++()
    {
        ++scalar_column_index;
        ++simd_column_index;
        return *this;
    }

    //!\brief Increments the iterator to the next column and returns the iterator pointing to the previous one.
    coordinate_matrix_simd_iterator operator++(int)
    {
        coordinate_matrix_simd_iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(coordinate_matrix_simd_iterator const & lhs, coordinate_matrix_simd_iterator const & rhs)
    {
        return lhs.scalar_column_index == rhs.scalar_column_index;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(coordinate_matrix_simd_iterator const & lhs, coordinate_matrix_simd_iterator const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}
};

} // namespace seqan3::detail
