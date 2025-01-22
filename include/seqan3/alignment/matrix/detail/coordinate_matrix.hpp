// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::coordinate_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/container/aligned_allocator.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/simd/views/iota_simd.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

namespace seqan3::detail
{

//------------------------------------------------------------------------------
// coordinate_matrix
//------------------------------------------------------------------------------

/*!\brief A matrix over coordinates.
 * \ingroup alignment_matrix
 * \implements std::ranges::forward_range
 *
 * \tparam index_t The underlying matrix index type; must mode std::integral or seqan3::simd::simd_index.
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
    requires (std::integral<index_t> || simd_index<index_t>)
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
            return matrix_index<index_t>{row_index_type{row_index}, column_index_type{column_index}};
        }
    };

    //!\brief Type alias for the scalar type defined by the seqan3::simd::simd_traits type.
    template <typename simd_index_t>
    using lazy_scalar_type_t = typename simd_traits<simd_index_t>::scalar_type;

    //!\brief The internal size type which depends on `index_t` being a simd vector or a scalar type.
    using size_type = lazy_conditional_t<simd_concept<index_t>, lazy<lazy_scalar_type_t, index_t>, index_t>;

    // The coordinate matrix iterator.
    class iterator;

    //!\brief The number of columns (corresponds to the size of one row).
    size_type column_count{};
    //!\brief The number od rows (corresponds to the size of one column).
    size_type row_count{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    coordinate_matrix() = default;                                      //!< Defaulted.
    coordinate_matrix(coordinate_matrix const &) = default;             //!< Defaulted.
    coordinate_matrix(coordinate_matrix &&) = default;                  //!< Defaulted.
    coordinate_matrix & operator=(coordinate_matrix const &) = default; //!< Defaulted.
    coordinate_matrix & operator=(coordinate_matrix &&) = default;      //!< Defaulted.
    ~coordinate_matrix() = default;                                     //!< Defaulted.

    /*!\brief Resets the coordinate matrix with the given end column index and end row index representing the new
     *        dimensions of the matrix.
     *\tparam column_index_t The size type for the column count; must model std:integral.
     *\tparam row_index_t The size type for the row count; must model std:integral.
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
    template <std::integral column_index_t, std::integral row_index_t>
    void resize(column_index_type<column_index_t> const column_count,
                row_index_type<row_index_t> const row_count) noexcept
    {
        this->column_count = column_count.get();
        this->row_count = row_count.get();
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column of the matrix.
    iterator begin() const noexcept
    {
        return iterator{size_type{}, row_count};
    }

    //!\brief Returns the iterator pointing to the end column of the matrix.
    iterator end() const noexcept
    {
        return iterator{column_count, row_count};
    }
    //!\}
};

//------------------------------------------------------------------------------
// iterator
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
    requires (std::integral<index_t> || simd_index<index_t>)
class coordinate_matrix<index_t>::iterator
{
private:
    //!\brief The iota view type which depends on the index type.
    using iota_view_t = lazy_conditional_t<simd_index<index_t>,
                                           lazy<iota_simd_view, index_t>,
                                           decltype(std::views::iota(size_type{}, size_type{}))>;
    //!\brief The currently represented column index.
    index_t column_id{0};
    //!\brief The number of rows (corresponds to the size of a column).
    size_type row_count{0};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::declval<iota_view_t>()
                                | std::views::transform(convert_to_matrix_coordinate{index_t{} /*column_id*/}));
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
    iterator() = default;                             //!< Defaulted.
    iterator(iterator const &) = default;             //!< Defaulted.
    iterator(iterator &&) = default;                  //!< Defaulted.
    iterator & operator=(iterator const &) = default; //!< Defaulted.
    iterator & operator=(iterator &&) = default;      //!< Defaulted.
    ~iterator() = default;                            //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the current column index and the row index marking the end
     *        of the rows (size of one column).
     *
     * \param[in] column_id \copybrief seqan3::detail::coordinate_matrix::iterator::column_id
     * \param[in] row_count \copybrief seqan3::detail::coordinate_matrix::iterator::row_count
     */
    explicit iterator(size_type column_id, size_type row_count) noexcept : row_count{std::move(row_count)}
    {
        if constexpr (simd_index<index_t>)
            this->column_id = simd::fill<index_t>(std::move(column_id));
        else
            this->column_id = std::move(column_id);
    }
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    auto operator*() const
    {
        if constexpr (simd_index<index_t>)
        {
            return views::iota_simd<index_t>(size_type{}, row_count)
                 | std::views::transform(convert_to_matrix_coordinate{column_id});
        }
        else
        {
            return std::views::iota(size_type{}, row_count)
                 | std::views::transform(convert_to_matrix_coordinate{column_id});
        }
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator to the next column.
    iterator & operator++()
    {
        // clang: pre-increment of a SIMD vector does not work
        if constexpr (simd_index<index_t>)
            column_id = column_id + simd::fill<index_t>(1);
        else
            ++column_id;

        return *this;
    }

    //!\brief Increments the iterator to the next column and returns the iterator pointing to the previous column.
    iterator operator++(int)
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(iterator const & lhs, iterator const & rhs)
    {
        if constexpr (simd_index<index_t>)
            return lhs.column_id[0] == rhs.column_id[0];
        else
            return lhs.column_id == rhs.column_id;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(iterator const & lhs, iterator const & rhs)
    {
        return !(lhs == rhs);
    }
    //!\}
};
} // namespace seqan3::detail
