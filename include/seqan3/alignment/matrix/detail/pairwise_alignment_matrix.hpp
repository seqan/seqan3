// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <tuple>

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/matrix/detail/trace_cell_proxy.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

template <tuple_like tuple_t>
//!\cond
    requires (std::tuple_size_v<tuple_t> == 2)
//!\endcond
class pairwise_alignment_cell_proxy : public tuple_t
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_cell_proxy()= default; //!< Defaulted.
    pairwise_alignment_cell_proxy(pairwise_alignment_cell_proxy const &) = default; //!< Defaulted.
    pairwise_alignment_cell_proxy(pairwise_alignment_cell_proxy &&) = default; //!< Defaulted.
    pairwise_alignment_cell_proxy & operator=(pairwise_alignment_cell_proxy const &) = default; //!< Defaulted.
    pairwise_alignment_cell_proxy & operator=(pairwise_alignment_cell_proxy &&) = default; //!< Defaulted.
    ~pairwise_alignment_cell_proxy() = default; //!< Defaulted.

    // Inherit the tuple constructors and assignment.
    using tuple_t::tuple_t;
    using tuple_t::operator=;

    template <typename score_proxy_tuple_t, typename trace_proxy_tuple_t>
    pairwise_alignment_cell_proxy & operator=(std::pair<affine_cell_proxy<score_proxy_tuple_t>, trace_cell_proxy<trace_proxy_tuple_t>> const & other)
    {
        using std::get;
        get<0>(*this) = get<0>(other);
        get<1>(*this) = get<1>(other);
        return *this;
    }

    template <typename score_proxy_tuple_t, typename trace_proxy_tuple_t>
    pairwise_alignment_cell_proxy & operator=(std::pair<affine_cell_proxy<score_proxy_tuple_t>, trace_cell_proxy<trace_proxy_tuple_t>> && other)
    {
        using std::get;
        get<0>(*this) = std::move(get<0>(other));
        get<1>(*this) = std::move(get<1>(other));
        return *this;
    }
    //!\}

    /*!\name Optimal score
     * \{
     */
    //!\brief Access the optimal score of the wrapped score matrix cell.
    decltype(auto) optimal_score() & { using std::get; return get<0>(get<0>(*this)); }
    //!\overload
    decltype(auto) optimal_score() const & { using std::get; return get<0>(get<0>(*this)); }
    //!\overload
    decltype(auto) optimal_score() && { using std::get; return get<0>(std::move(get<0>(*this))); }
    //!\overload
    decltype(auto) optimal_score() const &&
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using std::get;
        return static_cast<std::tuple_element_t<0, tuple_t> const &&>(get<0>(std::move(get<0>(*this))));
    }
    //!\}

    /*!\name Horizontal score
     * \{
     */
    //!\brief Access the horizontal score of the wrapped score matrix cell.
    decltype(auto) horizontal_score() & { using std::get; return get<1>(get<0>(*this)); }
    //!\overload
    decltype(auto) horizontal_score() const & { using std::get; return get<1>(get<0>(*this)); }
    //!\overload
    decltype(auto) horizontal_score() && { using std::get; return get<1>(std::move(get<0>(*this))); }
    //!\overload
     decltype(auto) horizontal_score() const &&
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using std::get;
        return static_cast<std::tuple_element_t<1, tuple_t> const &&>(get<1>(std::move(get<0>(*this))));
    }
    //!\}

    /*!\name Vertical score
     * \{
     */
    //!\brief Access the vertical score of the wrapped score matrix cell.
    decltype(auto) vertical_score() & { using std::get; return get<2>(get<0>(*this)); }
    //!\overload
    decltype(auto) vertical_score() const & { using std::get; return get<2>(get<0>(*this)); }
    //!\overload
    decltype(auto) vertical_score() && { using std::get; return get<2>(std::move(get<0>(*this))); }
    //!\overload
    decltype(auto) vertical_score() const &&
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using std::get;
        return static_cast<std::tuple_element_t<2, tuple_t> const &&>(get<2>(std::move(get<0>(*this))));
    }

    /*!\name Trace
     * \{
     */
    //!\brief Access the optimal score of the wrapped score matrix cell.
    decltype(auto) trace() & { using std::get; return get<0>(get<1>(*this)); }
    //!\overload
    decltype(auto) trace() const & { using std::get; return get<0>(get<1>(*this)); }
    //!\overload
    decltype(auto) trace() && { using std::get; return get<0>(std::move(get<1>(*this))); }
    //!\overload
    decltype(auto) trace() const &&
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using std::get;
        return static_cast<std::tuple_element_t<0, tuple_t> const &&>(get<0>(std::move(get<1>(*this))));
    }
    //!\}

    /*!\name Horizontal trace
     * \{
     */
    //!\brief Access the horizontal score of the wrapped score matrix cell.
    decltype(auto) horizontal_trace() & { using std::get; return get<1>(get<1>(*this)); }
    //!\overload
    decltype(auto) horizontal_trace() const & { using std::get; return get<1>(get<1>(*this)); }
    //!\overload
    decltype(auto) horizontal_trace() && { using std::get; return get<1>(std::move(get<1>(*this))); }
    //!\overload
     decltype(auto) horizontal_trace() const &&
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using std::get;
        return static_cast<std::tuple_element_t<1, tuple_t> const &&>(get<1>(std::move(get<1>(*this))));
    }
    //!\}

    /*!\name Vertical trace
     * \{
     */
    //!\brief Access the vertical score of the wrapped score matrix cell.
    decltype(auto) vertical_trace() & { using std::get; return get<2>(get<1>(*this)); }
    //!\overload
    decltype(auto) vertical_trace() const & { using std::get; return get<2>(get<1>(*this)); }
    //!\overload
    decltype(auto) vertical_trace() && { using std::get; return get<2>(std::move(get<1>(*this))); }
    //!\overload
    decltype(auto) vertical_trace() const &&
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using std::get;
        return static_cast<std::tuple_element_t<2, tuple_t> const &&>(get<2>(std::move(get<1>(*this))));
    }
    //!\}
};

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
template <typename score_matrix_t, typename trace_matrix_t>
class pairwise_alignment_matrix
{
private:

    using combined_matrix_t = decltype(views::zip(std::declval<score_matrix_t &>(), std::declval<trace_matrix_t &>()));

    using score_type = remove_cvref_t<std::tuple_element_t<0, range_innermost_value_t<score_matrix_t>>>;

    class iterator;
    class sentinel;

    score_matrix_t score_matrix{};
    trace_matrix_t trace_matrix{};
    combined_matrix_t combined_matrix{};

    //!\brief The number of columns for this matrix.
    size_t column_count{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_matrix() = default; //!< Defaulted.
    pairwise_alignment_matrix(pairwise_alignment_matrix const &) = default; //!< Defaulted.
    pairwise_alignment_matrix(pairwise_alignment_matrix &&) = default; //!< Defaulted.
    pairwise_alignment_matrix & operator=(pairwise_alignment_matrix const &) = default; //!< Defaulted.
    pairwise_alignment_matrix & operator=(pairwise_alignment_matrix &&) = default; //!< Defaulted.
    ~pairwise_alignment_matrix() = default; //!< Defaulted.
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
                row_index_type<row_index_t> const row_count,
                score_type const initial_score = score_type{})
    {
        score_matrix.resize(column_count, row_count, initial_score);
        trace_matrix.resize(column_count, row_count);
        combined_matrix = views::zip(score_matrix, trace_matrix);
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first column.
    iterator begin()
    {
        return iterator{combined_matrix.begin()};
    }

    //!\brief This score matrix is not const-iterable.
    iterator begin() const = delete;

    //!\brief Returns the sentinel pointing behind the last column.
    sentinel end()
    {
        return sentinel{combined_matrix.end()};
    }

    //!\brief This alignment matrix is not const-iterable.
    sentinel end() const = delete;
    //!\}

    /*!\brief Returns a trace path starting from the given coordinate and ending in the cell with
     *        seqan3::detail::trace_directions::none.
     * \param[in] trace_begin A seqan3::matrix_coordinate pointing to the begin of the trace to follow.
     * \returns A std::ranges::subrange over the corresponding trace path.
     * \throws std::invalid_argument if the specified coordinate is out of range.
     */
    auto trace_path(matrix_coordinate const & trace_begin)
    {
        return trace_matrix.trace_path(trace_begin);
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
template <typename score_matrix_t, typename trace_matrix_t>
class pairwise_alignment_matrix<score_matrix_t, trace_matrix_t>::iterator
{
private:
    //!\brief The type of the zipped score column.
    using combined_matrix_iterator_type = std::ranges::iterator_t<combined_matrix_t>;

    using score_matrix_reference_type = std::ranges::range_reference_t<score_matrix_t>;
    using trace_matrix_reference_type = std::ranges::range_reference_t<trace_matrix_t>;
    using combined_column_type = decltype(views::zip(std::declval<score_matrix_reference_type &&>(),
                                                     std::declval<trace_matrix_reference_type &&>()));

    //!\brief The transform adaptor to convert the tuple from the zip view into a seqan3::detail::affine_cell_type.
    static constexpr auto transform_to_combined_matrix_cell = std::views::transform([] (auto && tpl)
    {
        using fwd_tuple_t = decltype(tpl);
        return pairwise_alignment_cell_proxy<remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
    });

    //!\brief The pointer to the underlying matrix.
    combined_matrix_iterator_type combined_matrix_it{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = decltype(std::declval<combined_column_type>() | transform_to_combined_matrix_cell);
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
    explicit iterator(combined_matrix_iterator_type combined_matrix_it) noexcept :
        combined_matrix_it{std::move(combined_matrix_it)}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    combined_matrix_iterator_type base() const
        requires std::copy_constructible<combined_matrix_iterator_type>
    {
        return combined_matrix_it;
    }

    //!\brief Returns the range over the current column.
    reference operator*() const
    {
        auto && [score_column, trace_column] = *combined_matrix_it;
        return views::zip(std::move(score_column), std::move(trace_column)) | transform_to_combined_matrix_cell;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Move `this` to the next column.
    iterator & operator++()
    {
        ++combined_matrix_it;
        return *this;
    }

    //!\brief Move `this` to the next column.
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

};

template <typename score_matrix_t, typename trace_matrix_t>
class pairwise_alignment_matrix<score_matrix_t, trace_matrix_t>::sentinel
{
private:
    //!\brief The type of the zipped score column.
    using combined_matrix_sentinel_type = std::ranges::sentinel_t<combined_matrix_t>;

    //!\brief The pointer to the underlying matrix.
    combined_matrix_sentinel_type combined_matrix_sentinel{};
public:
    /*!\name Constructor, assignment and destructor
     * \{
     */
    sentinel() noexcept = default; //!< Defaulted.
    sentinel(sentinel const &) noexcept = default; //!< Defaulted.
    sentinel(sentinel &&) noexcept = default; //!< Defaulted.
    sentinel & operator=(sentinel const &) noexcept = default; //!< Defaulted.
    sentinel & operator=(sentinel &&) noexcept = default; //!< Defaulted.
    ~sentinel() = default; //!< Defaulted.

    /*!\brief Initialises the sentinel from the underlying matrix.
     *
     * \param[in] host_matrix The underlying matrix.
     * \param[in] initial_column_id The initial column index.
     */
    explicit sentinel(combined_matrix_sentinel_type combined_matrix_sentinel) noexcept :
        combined_matrix_sentinel{std::move(combined_matrix_sentinel)}
    {}
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(iterator const & lhs, sentinel const & rhs) noexcept
    {
        return lhs.base() == rhs.combined_matrix_sentinel;
    }

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(sentinel const & lhs, iterator const & rhs) noexcept
    {
        return rhs == lhs;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(iterator const & lhs, sentinel const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(sentinel const & lhs, iterator const & rhs) noexcept
    {
        return rhs != lhs;
    }
    //!\}
};

} // namespace seqan3::detail
