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

#include <seqan3/std/ranges>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>

namespace seqan3::detail
{

/*!\brief A function object that converts a column index and a row index to a seqan3::detail::matrix_coordinate.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This helper function object is used for the seqan3::detail::coordinate_matrix_iterator to convert a column index
 * and row index into a seqan3::detail::matrix_coordinate. The seqan3::detail::coordinate_matrix_iterator iterates
 * over the columns of the underlying coordinate matrix and returns a transformed iota range over the row indices of
 * this column. The column index is stored as member inside of this function object since every column index is
 * associated with many row indices.
 */
struct convert_to_matrix_coordinate
{
    //!\brief The index of the represented column.
    size_t column_index{0};

    /*!\brief The conversion operator
     *
     * \param[in] row_index The index of the represented row.
     *
     * \returns seqan3::detail::matrix_coordinate for the current column and row index.
     */
    auto operator()(size_t const row_index) noexcept
    {
        return matrix_coordinate{row_index_type{row_index}, column_index_type{column_index}};
    }
};

/*!\brief The iterator for the seqan3::detail::coordinate_matrix.
 * \ingroup alignment_matrix
 * \implements std::forward_iterator
 *
 * \details
 *
 * Iterates over the columns of the underlying seqan3::detail::coordinate_matrix. The iterator returns a transformed
 * std::ranges::views::iota view over the row indices. In the transformation every row index is converted to a
 * seqan3::detail::matrix_coordinate using the seqan3::detail::convert_to_matrix_coordinate function object.
 */
class coordinate_matrix_iterator
{
private:
    //!\brief The currently represented column index.
    size_t column_id{0};
    //!\brief The size of a column.
    size_t column_size{0};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::views::iota(size_t{0u}, size_t{1u})
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
     * \param[in] column_size \copybrief seqan3::detail::coordinate_matrix_iterator::column_size
     */
    explicit coordinate_matrix_iterator(size_t column_id, size_t column_size) noexcept :
        column_id{column_id},
        column_size{column_size}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    auto operator*() const
    {
        return std::views::iota(size_t{0u}, column_size)
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

/*!\brief A matrix over coordinates.
 * \ingroup alignment_matrix
 * \implements std::ranges::forward_range
 *
 * \details
 *
 * This matrix emulates a two-dimensional index, such that each cell inside of the alignment matrix can be located
 * through a unique coordinate (=index). In the alignment algorithm this matrix is paired with the alignment matrix
 * to form an indexed alignment matrix.
 *
 * This matrix is cheap as it stores only the dimensions of the matrix and does not allocate any memory for the
 * coordinates. It uses the seqan3::detail::matrix_coordinate_iterator to iterate
 * over the columns of this virtual matrix. When the seqan3::detail::matrix_coordinate_iterator is dereferenced it
 * returns an on-the-fly constructed range ([prvalue](https://en.cppreference.com/w/cpp/language/value_category))
 * representing the current indexed column. The reference type of this range is seqan3::detail::matrix_coordinate.
 */
class coordinate_matrix
{
private:
    //!\brief Represents the size of a row (number of columns).
    size_t number_of_columns{};
    //!\brief Represents the size of a column (number of rows).
    size_t number_of_rows{};

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

    /*!\brief Resets the coordinate matrix with the given end column index and end row index representing the dimensions
     *        of the reset matrix.
     *
     * \param[in] number_of_columns \copybrief seqan3::detail::coordinate_matrix::number_of_columns
     * \param[in] number_of_rows \copybrief seqan3::detail::coordinate_matrix::number_of_rows
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
    void resize(column_index_type<size_t> const number_of_columns, row_index_type<size_t> const number_of_rows) noexcept
    {
        this->number_of_columns = number_of_columns.get();
        this->number_of_rows = number_of_rows.get();
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column of the matrix.
    coordinate_matrix_iterator begin() const noexcept
    {
        return coordinate_matrix_iterator{0u, number_of_rows};
    }

    //!\brief Returns the iterator pointing to the end column of the matrix.
    coordinate_matrix_iterator end() const noexcept
    {
        return coordinate_matrix_iterator{number_of_columns, number_of_rows};
    }
    //!\}
};

} // namespace seqan3::detail
