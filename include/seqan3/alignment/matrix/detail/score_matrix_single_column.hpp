// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::score_matrix_single_column.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <vector>

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/container/aligned_allocator.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/views/repeat_n.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief Score matrix for the pairwise alignment using only a single column.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam score_t The type of the score; must model seqan3::arithmetic or seqan3::simd::simd_concept.
 *
 * \details
 *
 * In many cases it is sufficient to store only a single score column to compute the alignment between two sequences.
 * Since the alignment is computed iteratively column by column, the same memory can be reused for the next score.
 * This score matrix stores the complete column for both the optimal and horizontal score, but only stores a single
 * value for the vertical column. Hence, this matrix can only be used for a column
 * based computation layout.
 *
 * ### Range interface
 *
 * The matrix offers a input range interface over the columns of the matrix. Dereferencing the iterator will return
 * another range which represents the actual score column in memory. The returned range is a
 * transformed seqan3::views::zip view over the optimal, horizontal and vertical column. The reference type of this
 * view is the seqan3::detail::affine_cell_proxy, which offers a practical interface to access the value of the
 * optimal, horizontal and vertical value of the underlying matrices.
 */
template <typename score_t>
    requires (arithmetic<score_t> || simd_concept<score_t>)
class score_matrix_single_column
{
private:
    //!\brief The type of the score column which allocates memory for the entire column.
    using physical_column_t = std::vector<score_t, aligned_allocator<score_t, alignof(score_t)>>;
    //!\brief The type of the virtual score column which only stores one value.
    using virtual_column_t = decltype(views::repeat_n(score_t{}, 1));

    class matrix_iterator;

    //!\brief The column over the optimal scores.
    physical_column_t optimal_column{};
    //!\brief The column over the horizontal gap scores.
    physical_column_t horizontal_column{};
    //!\brief The virtual column over the vertical gap scores.
    virtual_column_t vertical_column{};
    //!\brief The number of columns for this matrix.
    size_t number_of_columns{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    score_matrix_single_column() = default;                                               //!< Defaulted.
    score_matrix_single_column(score_matrix_single_column const &) = default;             //!< Defaulted.
    score_matrix_single_column(score_matrix_single_column &&) = default;                  //!< Defaulted.
    score_matrix_single_column & operator=(score_matrix_single_column const &) = default; //!< Defaulted.
    score_matrix_single_column & operator=(score_matrix_single_column &&) = default;      //!< Defaulted.
    ~score_matrix_single_column() = default;                                              //!< Defaulted.

    //!\}

    /*!\brief Resizes the matrix.
     * \tparam column_index_t The column index type; must model std::integral.
     * \tparam row_index_t The row index type; must model std::integral.
     *
     * \param[in] number_of_columns The number of columns for this matrix.
     * \param[in] number_of_rows The number of rows for this matrix.
     * \param[in] initial_value Optional initial score value to use when resizing the underlying container.
     *
     * \details
     *
     * Resizes the optimal and the horizontal score column to the given number of rows and stores the number of
     * columns to created a counted iterator over the matrix columns.
     * Note the alignment matrix requires the number of columns and rows to be one bigger than the size of sequence1,
     * respectively sequence2.
     * Reallocation happens only if the new column size exceeds the current capacity of the optimal and horizontal
     * score column. The underlying vectors are initialised with the given `initial_value` or the default value of
     * the class's score type.
     *
     * ### Complexity
     *
     * Linear in the number of rows.
     *
     * ### Exception
     *
     * Basic exception guarantee. Might throw std::bad_alloc on resizing the internal columns.
     */
    template <std::integral column_index_t, std::integral row_index_t>
    void resize(column_index_type<column_index_t> const number_of_columns,
                row_index_type<row_index_t> const number_of_rows,
                score_t const initial_value = score_t{})
    {
        this->number_of_columns = number_of_columns.get();
        optimal_column.clear();
        horizontal_column.clear();
        SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(-Wstringop-overflow)
        optimal_column.resize(number_of_rows.get(), initial_value);
        SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_STOP
        horizontal_column.resize(number_of_rows.get(), initial_value);
        vertical_column = views::repeat_n(initial_value, number_of_rows.get());
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column.
    matrix_iterator begin()
    {
        return matrix_iterator{*this, 0u};
    }

    //!\brief This score matrix is not const-iterable.
    matrix_iterator begin() const = delete;

    //!\brief Returns the iterator pointing behind the last column.
    matrix_iterator end()
    {
        return matrix_iterator{*this, number_of_columns};
    }

    //!\brief This score matrix is not const-iterable.
    matrix_iterator end() const = delete;
    //!\}
};

/*!\brief Score matrix iterator for the pairwise alignment using only a single column.
 * \implements std::input_iterator
 *
 * \details
 *
 * Implements a counted iterator to simulate the iteration over the actual matrix. When dereferenced, the
 * iterator returns a view over the allocated memory of the respective columns. The returned view zips
 * the three columns into a single range and transforms the returned tuple to a
 * seqan3::detail::affine_cell_proxy to simplify the access to the correct values without knowing the internal
 * tuple layout returned by the seqan3::views::zip view.
 */
template <typename score_t>
    requires (arithmetic<score_t> || simd_concept<score_t>)
class score_matrix_single_column<score_t>::matrix_iterator
{
private:
    //!\brief The type of the zipped score column.
    using matrix_column_t = decltype(views::zip(std::declval<physical_column_t &>(),
                                                std::declval<physical_column_t &>(),
                                                std::declval<virtual_column_t &>()));

    //!\brief The transform adaptor to convert the tuple from the zip view into a seqan3::detail::affine_cell_type.
    static constexpr auto transform_to_affine_cell = std::views::transform(
        [](auto && tpl) -> affine_cell_proxy<std::remove_cvref_t<decltype(tpl)>>
        {
            using fwd_tuple_t = decltype(tpl);
            return affine_cell_proxy<std::remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
        });

    //!\brief The pointer to the underlying matrix.
    score_matrix_single_column * host_ptr{nullptr};
    //!\brief The current column index.
    size_t current_column_id{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::declval<matrix_column_t>() | transform_to_affine_cell);
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    matrix_iterator() noexcept = default;                                    //!< Defaulted.
    matrix_iterator(matrix_iterator const &) noexcept = default;             //!< Defaulted.
    matrix_iterator(matrix_iterator &&) noexcept = default;                  //!< Defaulted.
    matrix_iterator & operator=(matrix_iterator const &) noexcept = default; //!< Defaulted.
    matrix_iterator & operator=(matrix_iterator &&) noexcept = default;      //!< Defaulted.
    ~matrix_iterator() = default;                                            //!< Defaulted.

    /*!\brief Initialises the iterator from the underlying matrix.
     *
     * \param[in] host_matrix The underlying matrix.
     * \param[in] initial_column_id The initial column index.
     */
    explicit matrix_iterator(score_matrix_single_column & host_matrix, size_t const initial_column_id) noexcept :
        host_ptr{std::addressof(host_matrix)},
        current_column_id{initial_column_id}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the range over the current column.
    reference operator*() const
    {
        return views::zip(host_ptr->optimal_column, host_ptr->horizontal_column, host_ptr->vertical_column)
             | transform_to_affine_cell;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Move `this` to the next column.
    matrix_iterator & operator++()
    {
        ++current_column_id;
        return *this;
    }

    //!\brief Move `this` to the next column.
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(matrix_iterator const & lhs, matrix_iterator const & rhs) noexcept
    {
        return lhs.current_column_id == rhs.current_column_id;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(matrix_iterator const & lhs, matrix_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};

} // namespace seqan3::detail
