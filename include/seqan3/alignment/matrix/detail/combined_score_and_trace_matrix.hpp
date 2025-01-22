// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::combined_score_and_trace_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{
/*!\brief An alignment matrix that combines a score matrix with a trace matrix into a common interface.
 * \ingroup alignment_matrix
 * \implements std::ranges::input_range
 *
 * \tparam score_matrix_t The type of the score matrix; must model std::ranges::input_range (for more requirements see
 *                        details section).
 * \tparam trace_matrix_t The type of the trace matrix; must model std::ranges::input_range (for more requirements see
 *                        details section).
 *
 * \details
 *
 * To compute the actual, alignment an additional trace matrix is required. It stores the trace directions
 * to the previous cells for any computed cell within the alignment matrix. The following matrix combines a given score
 * with a trace matrix into a complete alignment matrix. The iterator increments over the columns of each underlying
 * matrix and combines them into a zip view, i.e. the combined alignment matrix is a range over ranges with the inner
 * range also modelling std::ranges::input_range.
 * To provide a uniform interface to the alignment algorithm, the returned zip view is further transformed such that
 * its reference type is a seqan3::detail::affine_cell_proxy.
 */
template <std::ranges::input_range score_matrix_t, std::ranges::input_range trace_matrix_t>
    requires (std::ranges::input_range<std::ranges::range_reference_t<score_matrix_t>>
              && std::ranges::input_range<std::ranges::range_reference_t<trace_matrix_t>>)
class combined_score_and_trace_matrix
{
private:
    //!\brief The score type extracted from the score matrix. (TODO: make this less complicated by exposing the score type).
    using score_type = std::remove_cvref_t<std::tuple_element_t<0, range_innermost_value_t<score_matrix_t>>>;

    class iterator;
    class sentinel;

    //!\brief The underlying score matrix.
    score_matrix_t score_matrix{};
    //!\brief The underlying trace matrix.
    trace_matrix_t trace_matrix{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    combined_score_and_trace_matrix() = default;                                                    //!< Defaulted.
    combined_score_and_trace_matrix(combined_score_and_trace_matrix const &) = default;             //!< Defaulted.
    combined_score_and_trace_matrix(combined_score_and_trace_matrix &&) = default;                  //!< Defaulted.
    combined_score_and_trace_matrix & operator=(combined_score_and_trace_matrix const &) = default; //!< Defaulted.
    combined_score_and_trace_matrix & operator=(combined_score_and_trace_matrix &&) = default;      //!< Defaulted.
    ~combined_score_and_trace_matrix() = default;                                                   //!< Defaulted.

    //!\}

    /*!\brief Resizes the matrix.
     * \tparam column_index_t The column index type; must model std::integral.
     * \tparam row_index_t The row index type; must model std::integral.
     *
     * \param[in] column_count The number of columns for this matrix.
     * \param[in] row_count The number of rows for this matrix.
     * \param[in] initial_score The initial score used to initialise the score matrix.
     *
     * \details
     *
     * Resizes the underlying score and trace matrix to the given dimensions.
     *
     * ### Complexity
     *
     * The complexity depends on the resize complexity of the underlying score and trace matrix.
     *
     * ### Exception
     *
     * Strong exception guarantee. Might throw std::bad_alloc.
     */
    template <std::integral column_index_t, std::integral row_index_t>
    void resize(column_index_type<column_index_t> const column_count,
                row_index_type<row_index_t> const row_count,
                score_type const initial_score = score_type{})
    {
        score_matrix_t tmp_score_matrix{};
        tmp_score_matrix.resize(column_count, row_count, initial_score);

        trace_matrix_t tmp_trace_matrix{};
        tmp_trace_matrix.resize(column_count, row_count);

        score_matrix = std::move(tmp_score_matrix);
        trace_matrix = std::move(tmp_trace_matrix);
    }

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator pointing to the first alignment column.
    iterator begin()
    {
        return iterator{score_matrix.begin(), trace_matrix.begin()};
    }

    //!\brief This score matrix is not const-iterable.
    iterator begin() const = delete;

    //!\brief Returns the sentinel pointing behind the last alignment column.
    sentinel end()
    {
        return sentinel{score_matrix.end()};
    }

    //!\brief This alignment matrix is not const-iterable.
    sentinel end() const = delete;
    //!\}

    /*!\brief Returns a trace path starting from the given coordinate and ending in the cell with
     *        seqan3::detail::trace_directions::none.
     * \param[in] from_coordinate A seqan3::matrix_coordinate pointing to the start of the trace to follow.
     *
     * \returns A std::ranges::subrange over the corresponding trace path.
     *
     * \throws std::invalid_argument if the specified coordinate is out of range.
     */
    auto trace_path(matrix_coordinate const & from_coordinate) const
    {
        return trace_matrix.trace_path(from_coordinate);
    }
};

/*!\brief Combined score and trace matrix iterator for the pairwise sequence alignment.
 * \implements std::input_iterator
 *
 * \details
 *
 * Implements a zipped iterator over the given score and trace matrix. It assumes that both matrices have the
 * same number of columns. If this invariant does not hold, using this iterator is undefined behaviour.
 * When dereferencing the iterator, the current score matrix column and trace matrix column will be zipped using the
 * seqan3::views::zip view and additionally transformed such that the reference type of the returned alignment column
 * is a seqan3::detail::affine_cell_proxy holding the score and trace values of the respective alignment matrix cells.
 */
template <std::ranges::input_range score_matrix_t, std::ranges::input_range trace_matrix_t>
    requires (std::ranges::input_range<std::ranges::range_reference_t<score_matrix_t>>
              && std::ranges::input_range<std::ranges::range_reference_t<trace_matrix_t>>)
class combined_score_and_trace_matrix<score_matrix_t, trace_matrix_t>::iterator
{
private:
    //!\brief The reference type of the score matrix.
    using score_matrix_reference_type = std::ranges::range_reference_t<score_matrix_t>;
    //!\brief The reference type of the trace matrix.
    using trace_matrix_reference_type = std::ranges::range_reference_t<trace_matrix_t>;

    static_assert(std::ranges::viewable_range<score_matrix_reference_type>);
    static_assert(std::ranges::viewable_range<trace_matrix_reference_type>);

    //!\brief The combined column type.
    using combined_column_type =
        decltype(views::zip(std::declval<score_matrix_reference_type>(), std::declval<trace_matrix_reference_type>()));
    //!\brief The type of the score matrix iterator.
    using score_matrix_iter_type = std::ranges::iterator_t<score_matrix_t>;
    //!\brief The type of the trace matrix iterator.
    using trace_matrix_iter_type = std::ranges::iterator_t<trace_matrix_t>;

    // Befriend the base class.
    template <std::ranges::input_range other_score_matrix_t, std::ranges::input_range other_trace_matrix_t>
        requires (std::ranges::input_range<std::ranges::range_reference_t<other_score_matrix_t>>
                  && std::ranges::input_range<std::ranges::range_reference_t<other_trace_matrix_t>>)
    friend class combined_score_and_trace_matrix;

    //!\brief The transform adaptor to convert the tuple from the zip view into a seqan3::detail::affine_cell_type.
    static constexpr auto transform_to_combined_matrix_cell = std::views::transform(
        [](auto && tpl) -> affine_cell_proxy<std::remove_cvref_t<decltype(tpl)>>
        {
            using fwd_tuple_t = decltype(tpl);
            return affine_cell_proxy<std::remove_cvref_t<fwd_tuple_t>>{std::forward<fwd_tuple_t>(tpl)};
        });

    //!\brief The current iterator over the score matrix.
    score_matrix_iter_type score_matrix_it;
    //!\brief The current iterator over the trace matrix.
    trace_matrix_iter_type trace_matrix_it;

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
    //!\brief The iterator concept.
    using iterator_concept = std::input_iterator_tag;
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

    /*!\brief Initialises the iterator from the underlying matrix.
     *
     * \param[in] score_matrix_it \copybrief seqan3::detail::combined_score_and_trace_matrix::iterator::score_matrix_it
     * \param[in] trace_matrix_it \copybrief seqan3::detail::combined_score_and_trace_matrix::iterator::trace_matrix_it
     */
    explicit iterator(score_matrix_iter_type score_matrix_it, trace_matrix_iter_type trace_matrix_it) noexcept :
        score_matrix_it{std::move(score_matrix_it)},
        trace_matrix_it{std::move(trace_matrix_it)}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the range over the current column.
    reference operator*() const
    {
        return views::zip(*score_matrix_it, *trace_matrix_it) | transform_to_combined_matrix_cell;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance the iterator to the next alignment matrix column.
    iterator & operator++()
    {
        ++score_matrix_it;
        ++trace_matrix_it;
        return *this;
    }

    //!\brief Advance the iterator to the next alignment matrix column.
    void operator++(int)
    {
        ++(*this);
    }
    //!\}
};

/*!\brief The sentinel type for the seqan3::detail::combined_score_and_trace_matrix.
 *
 * \details
 *
 * The sentinel only compares against the end of the score matrix and not the trace matrix,
 * since the invariant of the matrix requires that both submatrices have the same number of columns.
 */
template <std::ranges::input_range score_matrix_t, std::ranges::input_range trace_matrix_t>
    requires (std::ranges::input_range<std::ranges::range_reference_t<score_matrix_t>>
              && std::ranges::input_range<std::ranges::range_reference_t<trace_matrix_t>>)
class combined_score_and_trace_matrix<score_matrix_t, trace_matrix_t>::sentinel
{
private:
    //!\brief The sentinel type of the underlying score matrix.
    using score_matrix_sentinel_type = std::ranges::sentinel_t<score_matrix_t>;

    //!\brief The sentinel of the score matrix.
    score_matrix_sentinel_type score_matrix_sentinel{};

public:
    /*!\name Constructor, assignment and destructor
     * \{
     */
    sentinel() = default;                             //!< Defaulted.
    sentinel(sentinel const &) = default;             //!< Defaulted.
    sentinel(sentinel &&) = default;                  //!< Defaulted.
    sentinel & operator=(sentinel const &) = default; //!< Defaulted.
    sentinel & operator=(sentinel &&) = default;      //!< Defaulted.
    ~sentinel() = default;                            //!< Defaulted.

    /*!\brief Initialises the sentinel from the underlying matrix.
     *
     * \param[in] score_matrix_sentinel \copydoc seqan3::detail::combined_score_and_trace_matrix::sentinel::score_matrix_sentinel
     */
    explicit sentinel(score_matrix_sentinel_type score_matrix_sentinel) noexcept :
        score_matrix_sentinel{std::move(score_matrix_sentinel)}
    {}
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Checks if the iterator reached the end of the matrix.
    friend bool operator==(iterator const & lhs, sentinel const & rhs) noexcept
    {
        return rhs.equal(lhs);
    }

    //!\brief Checks if the iterator reached the end of the matrix.
    friend bool operator==(sentinel const & lhs, iterator const & rhs) noexcept
    {
        return rhs == lhs;
    }

    //!\brief Checks if the iterator did not reach the end of the matrix.
    friend bool operator!=(iterator const & lhs, sentinel const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Checks if the iterator did not reach the end of the matrix.
    friend bool operator!=(sentinel const & lhs, iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

private:
    /*!\brief Compares the stored score matrix sentinel with the given iterator.
     *
     * \param[in] iter The iterator to compare with.
     *
     * \details
     *
     * The additional member function is used to access the private member of the
     * seqan3::detail::combined_score_and_trace_matrix::iterator from the sentinel.
     */
    constexpr bool equal(iterator const & iter) const noexcept
    {
        return iter.score_matrix_it == score_matrix_sentinel;
    }
};

} // namespace seqan3::detail
