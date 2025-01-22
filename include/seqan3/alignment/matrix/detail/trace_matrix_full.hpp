// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::trace_matrix_full.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>
#include <span>
#include <vector>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/container/aligned_allocator.hpp>
#include <seqan3/utility/views/repeat_n.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief Trace matrix for the pairwise alignment using the full trace matrix.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam trace_t The type of the trace; must be the same as seqan3::detail::trace_directions.
 *
 * \details
 *
 * In the default trace back implementation we allocate the entire matrix using one byte per cell to store the
 * seqan3::detail::trace_directions.
 *
 * ### Range interface
 *
 * The matrix offers an input range interface over the columns of the matrix. Dereferencing the iterator will return
 * another range which represents the actual trace column in memory. The returned range is a
 * seqan3::views::zip view over the current column referencing the best trace, as well as the horizontal and vertical
 * trace column.
 */
template <typename trace_t>
    requires std::same_as<trace_t, trace_directions>
class trace_matrix_full
{
private:
    //!\brief The type to store the complete trace matrix.
    using matrix_t =
        two_dimensional_matrix<trace_t, aligned_allocator<trace_t, sizeof(trace_t)>, matrix_major_order::column>;
    //!\brief The type of the score column which allocates memory for the entire column.
    using physical_column_t = std::vector<trace_t>;
    //!\brief The type of the virtual score column which only stores one value.
    using virtual_column_t = decltype(views::repeat_n(trace_t{}, 1));

    class iterator;

    //!\brief The full trace matrix.
    matrix_t complete_matrix{};
    //!\brief The column over the horizontal traces.
    physical_column_t horizontal_column{};
    //!\brief The virtual column over the vertical traces.
    virtual_column_t vertical_column{};
    //!\brief The number of columns for this matrix.
    size_t column_count{};
    //!\brief The number of rows for this matrix.
    size_t row_count{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    trace_matrix_full() = default;                                      //!< Defaulted.
    trace_matrix_full(trace_matrix_full const &) = default;             //!< Defaulted.
    trace_matrix_full(trace_matrix_full &&) = default;                  //!< Defaulted.
    trace_matrix_full & operator=(trace_matrix_full const &) = default; //!< Defaulted.
    trace_matrix_full & operator=(trace_matrix_full &&) = default;      //!< Defaulted.
    ~trace_matrix_full() = default;                                     //!< Defaulted.

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
     * Resizes the entire trace matrix storing the best trace path and the horizontal trace column.
     * Note the trace matrix requires the number of columns and rows to be one bigger than the size of sequence1,
     * respectively sequence2 for the initialisation of the matrix.
     * Reallocation happens only if the new column size exceeds the current capacity of the underlying trace matrix.
     *
     * ### Complexity
     *
     * In worst case `column_count` times `row_count` memory is allocated.
     *
     * ### Exception
     *
     * Basic exception guarantee. Might throw std::bad_alloc on resizing the internal matrices.
     */
    template <std::integral column_index_t, std::integral row_index_t>
    void resize(column_index_type<column_index_t> const column_count, row_index_type<row_index_t> const row_count)
    {
        this->column_count = column_count.get();
        this->row_count = row_count.get();
        complete_matrix.resize(number_rows{this->row_count}, number_cols{this->column_count});
        horizontal_column.resize(this->row_count);
        vertical_column = views::repeat_n(trace_t{}, this->row_count);
    }

    /*!\brief Returns a trace path starting from the given coordinate and ending in the cell with
     *        seqan3::detail::trace_directions::none.
     * \param[in] trace_begin A seqan3::matrix_coordinate pointing to the begin of the trace to follow.
     * \returns A std::ranges::subrange over the corresponding trace path.
     * \throws std::invalid_argument if the specified coordinate is out of range.
     */
    auto trace_path(matrix_coordinate const & trace_begin) const
    {
        using matrix_iter_t = std::ranges::iterator_t<matrix_t const>;
        using trace_iterator_t = trace_iterator<matrix_iter_t>;
        using path_t = std::ranges::subrange<trace_iterator_t, std::default_sentinel_t>;

        if (trace_begin.row >= row_count || trace_begin.col >= column_count)
            throw std::invalid_argument{"The given coordinate exceeds the matrix in vertical or horizontal direction."};

        return path_t{trace_iterator_t{complete_matrix.begin() + matrix_offset{trace_begin}}, std::default_sentinel};
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
};

/*!\brief Trace matrix iterator for the pairwise alignment using the full trace matrix.
 * \implements std::input_iterator
 *
 * \details
 *
 * Implements a counted iterator to keep track of the current column within the matrix. When dereferenced, the
 * iterator returns a view over the allocated memory of the respective columns. The returned view zips
 * the three columns into a single range.
 */
template <typename trace_t>
    requires std::same_as<trace_t, trace_directions>
class trace_matrix_full<trace_t>::iterator
{
private:
    //!\brief A lightweight representation of a single column from the complete matrix.
    using single_trace_column_type = std::span<trace_t>;
    //!\brief The type of the zipped score column.
    using matrix_column_type = decltype(views::zip(std::declval<single_trace_column_type>(),
                                                   std::declval<physical_column_t &>(),
                                                   std::declval<virtual_column_t &>()));
    //!\brief The column type as value type.
    using matrix_column_value_t = std::vector<std::ranges::range_value_t<matrix_column_type>>;

    // Defines a proxy that can be converted to the value type.
    class column_proxy;

    //!\brief The pointer to the underlying matrix.
    trace_matrix_full * host_ptr{nullptr};
    //!\brief The current column index.
    size_t current_column_id{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = matrix_column_value_t;
    //!\brief The reference type.
    using reference = column_proxy;
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
    iterator() noexcept = default;                             //!< Defaulted.
    iterator(iterator const &) noexcept = default;             //!< Defaulted.
    iterator(iterator &&) noexcept = default;                  //!< Defaulted.
    iterator & operator=(iterator const &) noexcept = default; //!< Defaulted.
    iterator & operator=(iterator &&) noexcept = default;      //!< Defaulted.
    ~iterator() = default;                                     //!< Defaulted.

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
        auto column_begin = host_ptr->complete_matrix.data() + current_column_id * host_ptr->row_count;
        single_trace_column_type single_trace_column{column_begin, column_begin + host_ptr->row_count};

        return column_proxy{
            views::zip(std::move(single_trace_column), host_ptr->horizontal_column, host_ptr->vertical_column)};
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

/*!\brief The proxy returned as reference type.
 * \implements std::ranges::view
 *
 * \details
 *
 * The proxy stores the column view of the current iterator and offers a dedicated conversion operator to
 * assign it to the value type of the iterator.
 */
template <typename trace_t>
    requires std::same_as<trace_t, trace_directions>
class trace_matrix_full<trace_t>::iterator::column_proxy : public std::ranges::view_interface<column_proxy>
{
private:
    //!\brief The represented column.
    matrix_column_type column;

public:
    /*!\name Constructor, assignment and destructor
     * \{
     */
    column_proxy() = default;                                 //!< Defaulted.
    column_proxy(column_proxy const &) = default;             //!< Defaulted.
    column_proxy(column_proxy &&) = default;                  //!< Defaulted.
    column_proxy & operator=(column_proxy const &) = default; //!< Defaulted.
    column_proxy & operator=(column_proxy &&) = default;      //!< Defaulted.
    ~column_proxy() = default;                                //!< Defaulted.

    /*!\brief Initialises the proxy with the respective column.
    *
    * \param[in] column The column to set.
    */
    explicit column_proxy(matrix_column_type && column) noexcept : column{std::move(column)}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns an iterator to the begin of the column.
    std::ranges::iterator_t<matrix_column_type> begin()
    {
        return column.begin();
    }
    //!\brief Const iterator is not accessible.
    std::ranges::iterator_t<matrix_column_type> begin() const = delete;

    //!\brief Returns a sentinel marking the end of the column.
    std::ranges::sentinel_t<matrix_column_type> end()
    {
        return column.end();
    }

    //!\brief Const sentinel is not accessible.
    std::ranges::sentinel_t<matrix_column_type> end() const = delete;
    //!\}

    //!\brief Implicitly converts the column proxy into the value type of the iterator.
    constexpr operator matrix_column_value_t() const
    {
        matrix_column_value_t target{};
        std::ranges::copy(column, std::back_inserter(target));
        return target;
    }
};

} // namespace seqan3::detail
