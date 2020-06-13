// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_matrix_full.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <seqan3/std/span>
#include <vector>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/trace_cell_proxy.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief Score matrix for the pairwise alignment using only a single column.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam score_t The type of the score; must model seqan3::arithmetic.
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
template <typename trace_t>
//!\cond
    requires std::same_as<trace_t, trace_directions>
//!\endcond
class trace_matrix_full
{
private:
    //!\brief The type to store the complete matrix.
    using matrix_t = two_dimensional_matrix<trace_t,
                                            aligned_allocator<trace_t, sizeof(trace_t)>,
                                            matrix_major_order::column>;
    //!\brief The type of the score column which allocates memory for the entire column.
    using physical_column_t = std::vector<trace_t>;
    //!\brief The type of the virtual score column which only stores one value.
    using virtual_column_t = decltype(views::repeat_n(trace_t{}, 1));

    class iterator;

    //!\brief The column over the optimal scores.
    matrix_t complete_matrix{};
    //!\brief The column over the horizontal gap scores.
    physical_column_t horizontal_column{};
    //!\brief The virtual column over the vertical gap scores.
    virtual_column_t vertical_column{};
    //!\brief The number of columns for this matrix.
    size_t column_count{};
    //!\brief The number of rows for this matrix.
    size_t row_count{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    trace_matrix_full() = default; //!< Defaulted.
    trace_matrix_full(trace_matrix_full const &) = default; //!< Defaulted.
    trace_matrix_full(trace_matrix_full &&) = default; //!< Defaulted.
    trace_matrix_full & operator=(trace_matrix_full const &) = default; //!< Defaulted.
    trace_matrix_full & operator=(trace_matrix_full &&) = default; //!< Defaulted.
    ~trace_matrix_full() = default; //!< Defaulted.
    //!\}

    /*!\brief Resizes the matrix.
     * \tparam column_index_t The column index type; must model std::integral.
     * \tparam row_index_t The row index type; must model std::integral.
     *
     * \param[in] column_count The number of columns for this matrix.
     * \param[in] row_count The number of rows for this matrix.
     *
     * \details
     *
     * Resizes the optimal and the horizontal score column to the given number of rows and stores the number of
     * columns to created a counted iterator over the matrix columns.
     * Note the alignment matrix requires the number of columns and rows to be one bigger than the size of sequence1,
     * respectively sequence2.
     * Reallocation happens only if the new column size exceeds the current capacity of the optimal and horizontal
     * score column. The underlying vectors will not be cleared during the reset.
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
    void resize(column_index_type<column_index_t> const column_count,
                row_index_type<row_index_t> const row_count)
    {
        this->column_count = column_count.get();
        this->row_count = row_count.get();
        complete_matrix.resize(number_rows{this->row_count}, number_cols{this->column_count});
        horizontal_column.resize(this->row_count);
        vertical_column = views::repeat_n(trace_t{}, this->row_count);
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column.
    iterator begin()
    {
        return iterator{*this, 0u};
    }

    //!\brief This score matrix is not const-iterable.
    iterator begin() const = delete;

    //!\brief Returns the iterator pointing behind the last column.
    iterator end()
    {
        return iterator{*this, column_count};
    }

    //!\brief This score matrix is not const-iterable.
    iterator end() const = delete;
    //!\}

    /*!\brief Returns a trace path starting from the given coordinate and ending in the cell with
     *        seqan3::detail::trace_directions::none.
     * \param[in] trace_begin A seqan3::matrix_coordinate pointing to the begin of the trace to follow.
     * \returns A std::ranges::subrange over the corresponding trace path.
     * \throws std::invalid_argument if the specified coordinate is out of range.
     */
    auto trace_path(matrix_coordinate const & trace_begin)
    {
        using matrix_iter_t = std::ranges::iterator_t<matrix_t>;
        using trace_iterator_t = trace_iterator<matrix_iter_t>;
        using path_t = std::ranges::subrange<trace_iterator_t, std::default_sentinel_t>;

        if (trace_begin.row >= row_count || trace_begin.col >= column_count)
            throw std::invalid_argument{"The given coordinate exceeds the matrix in vertical or horizontal direction."};

        return path_t{trace_iterator_t{complete_matrix.begin() + matrix_offset{trace_begin}},
                      std::default_sentinel};
    }

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
template <typename trace_t>
//!\cond
    requires std::same_as<trace_t, trace_directions>
//!\endcond
class trace_matrix_full<trace_t>::iterator
{
private:
    //!\brief A lightweight representation of a single column from the complete matrix.
    using single_trace_column_t = std::span<trace_t>;
    //!\brief The type of the zipped score column.
    using matrix_column_t = decltype(views::zip(std::declval<single_trace_column_t>(),
                                                std::declval<physical_column_t &>(),
                                                std::declval<virtual_column_t &>()));

    //!\brief The transform adaptor to convert the tuple from the zip view into a seqan3::detail::trace_cell_type.
    static constexpr auto transform_to_trace_cell = std::views::transform([] (auto && tpl)
    {
        using fwd_tuple_t = decltype(tpl);
        return trace_cell_proxy<remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
    });

    //!\brief The pointer to the underlying matrix.
    trace_matrix_full * host_ptr{nullptr};
    //!\brief The current column index.
    size_t current_column_id{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::declval<matrix_column_t>() | transform_to_trace_cell);
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
    iterator() noexcept = default; //!< Defaulted.
    iterator(iterator const &) noexcept = default; //!< Defaulted.
    iterator(iterator &&) noexcept = default; //!< Defaulted.
    iterator & operator=(iterator const &) noexcept = default; //!< Defaulted.
    iterator & operator=(iterator &&) noexcept = default; //!< Defaulted.
    ~iterator() = default; //!< Defaulted.

    /*!\brief Initialises the iterator from the underlying matrix.
     *
     * \param[in] host_matrix The underlying matrix.
     * \param[in] initial_column_id The initial column index.
     */
    explicit iterator(trace_matrix_full & host_matrix, size_t const initial_column_id) noexcept :
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
        size_t const column_offset = current_column_id * host_ptr->row_count;
        single_trace_column_t single_trace_column{host_ptr->complete_matrix.data() + column_offset,
                                                  host_ptr->complete_matrix.data() + column_offset + host_ptr->row_count};
        return views::zip(std::move(single_trace_column), host_ptr->horizontal_column, host_ptr->vertical_column)
             | transform_to_trace_cell;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Move `this` to the next column.
    iterator & operator++()
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
    friend bool operator==(iterator const & lhs, iterator const & rhs) noexcept
    {
        return lhs.current_column_id == rhs.current_column_id;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(iterator const & lhs, iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};

} // namespace seqan3::detail
