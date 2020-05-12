// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::two_dimensional_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <memory>
#include <vector>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_base.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

//!\brief Strong type for setting the column dimension of a matrix.
//!\ingroup alignment_matrix
struct number_cols : strong_type<size_t, number_cols>
{
    //!\brief Import the base constructors.
    using strong_type<size_t, number_cols>::strong_type;
};

//!\brief Strong type for setting the row dimension of a matrix.
//!\ingroup alignment_matrix
struct number_rows : strong_type<size_t, number_rows>
{
    //!\brief Import the base constructors.
    using strong_type<size_t, number_rows>::strong_type;
};

/*!\brief A two dimensional matrix used inside of alignment algorithms.
 * \ingroup alignment_matrix
 * \implements std::ranges::random_access_range
 *
 * \tparam value_t The value type to store.
 * \tparam allocator_t The allocator type used to allocate the storage; defaults to std::allocator.
 * \tparam order The matrix-major-order; defaults to seqan3::detail::matrix_major_order::row.
 *
 * \details
 *
 * This two dimensional matrix type is a base data structure for several alignment matrices. It can be customised over
 * the value type, the  allocator type and the major matrix order. The data is stored in a flat std::vector.
 * Depending on the given `order`, the element access on the underlying buffer follows a row-major-order or a
 * column-major-order. Accordingly, it is cache friendly and thus more efficient to access the data in the
 * row-major-order row by row instead of column by column and vice versa for the column-major-order.
 */
template <typename value_t,
          typename allocator_t = std::allocator<value_t>,
          matrix_major_order order = matrix_major_order::row>
class two_dimensional_matrix
{
private:

    /*!\name Auxiliary types
     * \{
     */
    using storage_type = std::vector<value_t, allocator_t>; //!< The storage type.
    //!\}

    // Forward declaration. For definition see below.
    template <typename matrix_t>
    class iterator_type;

public:

    /*!\name Associated types
     * \{
     */
    using value_type = typename storage_type::value_type;  //!< The value type.
    using reference = typename storage_type::reference;  //!< The reference type.
    using const_reference = typename storage_type::const_reference; //!< The const reference type.
    using pointer = typename storage_type::pointer;  //!< The pointer type.
    using const_pointer = typename storage_type::const_pointer;  //!< The pointer type.
    using difference_type = typename storage_type::difference_type; //!< The difference type.
    using size_type = typename storage_type::size_type;  //!< The difference type.
    using iterator = iterator_type<two_dimensional_matrix>; //!< The iterator type.
    using const_iterator = iterator_type<two_dimensional_matrix const>; //!< The const iterator type.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    two_dimensional_matrix() = default; //!< Defaulted
    two_dimensional_matrix(two_dimensional_matrix const &) = default; //!< Defaulted
    two_dimensional_matrix(two_dimensional_matrix &&) = default; //!< Defaulted
    two_dimensional_matrix & operator=(two_dimensional_matrix const &) = default; //!< Defaulted
    two_dimensional_matrix & operator=(two_dimensional_matrix &&) = default; //!< Defaulted
    ~two_dimensional_matrix() = default; //!< Defaulted

    /*!\brief Constructs the matrix by the given dimensions.
     * \param row_dim The row dimension (number of rows).
     * \param col_dim The column dimension (number of columns).
     */
    two_dimensional_matrix(number_rows const row_dim, number_cols const col_dim) :
        row_dim{row_dim.get()}, col_dim{col_dim.get()}
    {
        storage.resize(row_dim.get() * col_dim.get());
    }

    /*!\brief Constructs the matrix by the given dimensions and initialises it with the given range.
     * \tparam entries_t Range of values that are convertible to value_type; must model std::ranges::forward_range
     * \param[in] row_dim The row dimension (number of rows).
     * \param[in] col_dim The column dimension (number of columns).
     * \param[in] entries A range used to fill the underlying matrix.
     */
    template <std::ranges::forward_range entries_t>
    //!\cond
        requires (std::convertible_to<std::ranges::range_value_t<entries_t>, value_type>)
    //!\endcond
    two_dimensional_matrix(number_rows const row_dim, number_cols const col_dim, entries_t entries) :
        row_dim{row_dim.get()},
        col_dim{col_dim.get()}
    {
        static_assert(std::move_constructible<std::ranges::range_value_t<entries_t>>, "The value type must be moveable.");

        assert(static_cast<size_t>(std::ranges::distance(entries)) == (row_dim.get() * col_dim.get()));
        storage.resize(row_dim.get() * col_dim.get());
        std::ranges::move(entries, storage.begin());
    }

    //!\overload
    two_dimensional_matrix(number_rows const row_dim, number_cols const col_dim, storage_type entries) :
        row_dim{row_dim.get()},
        col_dim{col_dim.get()}
    {
        assert(static_cast<size_t>(std::ranges::distance(entries)) == (row_dim.get() * col_dim.get()));
        storage = std::move(entries);
    }

    /*!\brief Explicit construction from the other major-order.
     * \tparam other_value_t The target value type; must be assignable from `value_t`.
     * \tparam other_allocator_t The allocator type used for the target matrix.
     * \tparam other_order The other seqan3::detail::matrix_major_order.
     *
     * \details
     *
     * Copies the matrix cell by cell, rearranging the stored elements in the internal memory to represent the
     * converted major-order.
     *
     * Consider the following matrix:
     *
     * 0  1  2  3
     * 4  5  6  7
     * 8  9  10 11
     *
     * In row-major-order the data is stored in a flat vector in the following way:
     * 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11
     *
     * Converting it to column-major-order will rearrange the elements:
     * 0, 4, 8, 1, 5, 9, 2, 6, 10, 3, 7, 11
     *
     * Note that the matrix is not transposed, so that the general layout as displayed above will remain the same.
     * It only changes the matrix major order, i.e. data stored row wise is now stored column wise and vice versa.
     */
    template <typename other_value_t, typename other_allocator_t, matrix_major_order other_order>
    //!\cond
        requires std::assignable_from<other_value_t &, value_t &>
    //!\endcond
    explicit constexpr two_dimensional_matrix(two_dimensional_matrix<other_value_t,
                                                                     other_allocator_t,
                                                                     other_order> const & matrix) :
        two_dimensional_matrix{number_rows{matrix.rows()}, number_cols{matrix.cols()}}
    {
        for (size_t i = 0; i < cols(); ++i)
        {
            for (size_t j = 0; j < rows(); ++j)
            {
                matrix_coordinate coord{row_index_type{j}, column_index_type{i}};
                (*this)[coord] = matrix[coord];
            }
        }
    }
    //!\}

    /*!\brief Returns a reference to the element at the given coordinate.
     * \param[in] coordinate The two-dimensional coordinate to access.
     */
    constexpr reference operator[](matrix_coordinate const & coordinate) noexcept
    {
        assert(coordinate.col >= 0u);
        assert(coordinate.row >= 0u);
        assert(coordinate.col < cols());
        assert(coordinate.row < rows());

        return *(begin() + matrix_offset{row_index_type{static_cast<std::ptrdiff_t>(coordinate.row)},
                                         column_index_type{static_cast<std::ptrdiff_t>(coordinate.col)}});
    }

    //!\copydoc operator[]
    constexpr const_reference operator[](matrix_coordinate const & coordinate) const noexcept
    {
        assert(coordinate.col >= 0);
        assert(coordinate.row >= 0);
        assert(static_cast<size_t>(coordinate.col) < cols());
        assert(static_cast<size_t>(coordinate.row) < rows());

        return *(begin() + matrix_offset{row_index_type{static_cast<std::ptrdiff_t>(coordinate.row)},
                                         column_index_type{static_cast<std::ptrdiff_t>(coordinate.col)}});
    }

    //!\copydoc seqan3::detail::matrix::at
    constexpr reference at(matrix_coordinate const & coordinate)
    {
        if (static_cast<size_t>(coordinate.col) >= cols())
            throw std::invalid_argument{"Column index is out of range. Please check the dimensions of the matrix."};
        if (static_cast<size_t>(coordinate.row) >= rows())
            throw std::invalid_argument{"Row index is out of range. Please check the dimensions of the matrix."};

        return (*this)[coordinate];
    }

    //!\copydoc seqan3::detail::matrix::at
    constexpr const_reference at(matrix_coordinate const & coordinate) const
    {
        if (static_cast<size_t>(coordinate.col) >= cols())
            throw std::invalid_argument{"Column index is out of range. Please check the dimensions of the matrix."};
        if (static_cast<size_t>(coordinate.row) >= rows())
            throw std::invalid_argument{"Row index is out of range. Please check the dimensions of the matrix."};

        return (*this)[coordinate];
    }

    //!\copydoc seqan3::detail::matrix::rows
    size_t rows() const noexcept
    {
        return row_dim;
    }

    //!\copydoc seqan3::detail::matrix::cols
    size_t cols() const noexcept
    {
        return col_dim;
    }

    //!\brief Returns a pointer to the data.
    constexpr pointer data() noexcept
    {
        return storage.data();
    }

    //!\brief Returns a pointer to the data.
    constexpr const_pointer data() const noexcept
    {
        return storage.data();
    }

    /*!\name Iterators
     * \brief Provides iterators to move in two-dimensional space.
     * \{
     */
    //!\brief Returns an iterator pointing to the first element of the matrix.
    constexpr iterator begin() noexcept
    {
        return {*this, storage.begin()};
    }
    //!\copydoc two_dimensional_matrix::begin
    constexpr const_iterator begin() const noexcept
    {
        return {*this, storage.begin()};
    }

    //!\copydoc two_dimensional_matrix::begin
    constexpr const_iterator cbegin() const noexcept
    {
        return begin();
    }

    //!\brief Returns an iterator pointing behind-the-end of the matrix.
    constexpr iterator end() noexcept
    {
        return {*this, storage.end()};
    }

    //!\copydoc two_dimensional_matrix::end
    constexpr const_iterator end() const noexcept
    {
        return {*this, storage.end()};
    }

    //!\copydoc two_dimensional_matrix::end
    constexpr const_iterator cend() const noexcept
    {
        return end();
    }
    //!\}

private:

    storage_type storage; //!< The matrix as a one-dimensional (flattened) vector of entries.
    size_type    row_dim; //!< \copydoc seqan3::detail::matrix::rows
    size_type    col_dim; //!< \copydoc seqan3::detail::matrix::cols
};

/*!\brief A two-dimensional matrix iterator.
 * \tparam matrix_t The underlying seqan3::detail::two_dimensional_matrix.
 *
 * \details
 *
 * Offers a two-dimensional iterator interface over the seqan3::detail::two_dimensional_matrix, which stores the
 * data in a flattened one-dimensional vector.
 */
template <typename value_t, typename allocator_t, matrix_major_order order>
template <typename matrix_t>
class two_dimensional_matrix<value_t, allocator_t, order>::iterator_type :
    public two_dimensional_matrix_iterator_base<iterator_type<matrix_t>, order>
{
private:
    //!\brief The base class type.
    using base_t = two_dimensional_matrix_iterator_base<iterator_type<matrix_t>, order>;

    //!\brief Befriend the base crtp class.
    template <typename derived_t, matrix_major_order other_order>
    //!\cond
        requires is_type_specialisation_of_v<derived_t, iterator_type> && other_order == order
    //!\endcond
    friend class two_dimensional_matrix_iterator_base;

    //!\brief Befriend the corresponding const iterator.
    template <typename other_matrix_t>
        requires std::same_as<other_matrix_t, std::remove_const_t<matrix_t>> &&
                    std::is_const_v<matrix_t>
    friend class iterator_type;

    //!\brief The iterator of the underlying storage.
    using storage_iterator = std::ranges::iterator_t<
                                std::conditional_t<
                                    std::is_const_v<matrix_t>, storage_type const, storage_type
                                >
                            >;

public:

    /*!\name Associated types
    * \{
    */
    //!\brief Value type of this iterator.
    using value_type = std::iter_value_t<storage_iterator>;
    //!\brief Reference to `value_type`.
    using reference = std::iter_reference_t<storage_iterator>;
    //!\brief The pointer type.
    using pointer = typename storage_iterator::pointer;
    //!\brief Type for distances between iterators.
    using difference_type = std::iter_difference_t<storage_iterator>;
    //!\brief The iterator tag.
    using iterator_category = std::random_access_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
    * \{
    */
    constexpr iterator_type() = default; //!< Defaulted.
    constexpr iterator_type(iterator_type const &) = default; //!< Defaulted.
    constexpr iterator_type(iterator_type &&) = default; //!< Defaulted.
    constexpr iterator_type & operator=(iterator_type const &) = default; //!< Defaulted.
    constexpr iterator_type & operator=(iterator_type &&) = default; //!< Defaulted.
    ~iterator_type() = default; //!< Defaulted.

    /*!\brief Construction from the underlying matrix and the iterator over actual storage.
    * \param[in] matrix The underlying matrix to access the corresponding dimensions.
    * \param[in] iter   The underlying iterator over the actual storage; must model std::random_access_iterator.
    */
    constexpr iterator_type(matrix_t & matrix, storage_iterator iter) :
        matrix_ptr{&matrix},
        host_iter{iter}
    {}

    //!\brief Construction of cons-iterator from non-const-iterator.
    template <typename other_matrix_t>
    //!\cond
        requires std::same_as<other_matrix_t, std::remove_const_t<matrix_t>> && std::is_const_v<matrix_t>
    //!\endcond
    constexpr iterator_type(
        iterator_type<other_matrix_t> other) noexcept :
        matrix_ptr{other.matrix_ptr},
        host_iter{other.host_iter}
    {}
    //!\}

    // Import advance operator from base class.
    using base_t::operator+=;

    //!\brief Advances the iterator by the given `offset`.
    constexpr iterator_type & operator+=(matrix_offset const & offset) noexcept
    {
        assert(matrix_ptr != nullptr);

        if constexpr (order == matrix_major_order::column)
        {
            host_iter += (offset.col * matrix_ptr->rows());
            host_iter += offset.row;
        }
        else
        {
            host_iter += offset.col;
            host_iter += (offset.row * matrix_ptr->cols());
        }
        return *this;
    }

    //!\copydoc seqan3::detail::two_dimensional_matrix_iterator::coordinate()
    matrix_coordinate coordinate() const noexcept
    {
        assert(matrix_ptr != nullptr);

        auto diff = *this - matrix_ptr->begin();
        if constexpr (order == matrix_major_order::column)
            return {row_index_type{diff % matrix_ptr->rows()}, column_index_type{diff / matrix_ptr->rows()}};
        else
            return {row_index_type{diff / matrix_ptr->cols()}, column_index_type{diff % matrix_ptr->cols()}};
    }

private:

    matrix_t * matrix_ptr{nullptr};  //!< Points to the associated matrix.
    storage_iterator host_iter{}; //!< The underlying storage iterator.
};

} // namespace seqan3::detail
