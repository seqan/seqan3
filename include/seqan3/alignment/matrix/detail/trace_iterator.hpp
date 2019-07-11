// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_concept.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_base.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A trace iterator for affine gaps.
 * \ingroup alignment_matrix
 * \extends std::ForwardIterator
 *
 * \tparam matrix_iter_t The wrapped matrix iterator; must model seqan3::detail::TwoDimensionalMatrixIterator and
 *                       the iterator must be a trace matrix, i.e.
 *                       `std::Same<value_type_t<matrix_iter_t>, trace_directions>` must evaluate to `true`.
 *
 * \details
 *
 * This iterator follows the trace path as computed for affine gaps.
 * When dereferencing it outputs the solely direction seqan3::detail::trace_directions::diagonal,
 * seqan3::detail::trace_directions::up, or seqan3::detail::trace_directions::left. It does not directly
 * dereference the actual trace direction stored in the underlying matrix. Thus, it cannot be used as an output
 * iterator. When advancing this iterator it actually moves from right to left and from bottom to top in the
 * underlying matrix until an entry with seqan3::detail::trace_directions::none is found.
 */
template <TwoDimensionalMatrixIterator matrix_iter_t>
class trace_iterator
{
private:

    static_assert(std::Same<value_type_t<matrix_iter_t>, trace_directions>,
                  "Value type of the underlying iterator must be trace_directions");

    //!\brief Befriend with iterator over const verion.
    template <TwoDimensionalMatrixIterator other_matrix_iter_t>
        requires std::Constructible<matrix_iter_t, other_matrix_iter_t>
    friend class trace_iterator;

public:
    /*!\name Associated types
     * \{
     */
    using value_type = trace_directions; //!< The value type.
    using reference = value_type; //!< The reference type.
    using pointer = void; //!< The pointer type.
    using difference_type = std::ptrdiff_t; //!< The difference type.
    using iterator_category = std::forward_iterator_tag; //!< Forward iterator tag.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr trace_iterator() = default; //!< Defaulted.
    constexpr trace_iterator(trace_iterator const &) = default; //!< Defaulted.
    constexpr trace_iterator(trace_iterator &&) = default; //!< Defaulted.
    constexpr trace_iterator & operator=(trace_iterator const &) = default; //!< Defaulted.
    constexpr trace_iterator & operator=(trace_iterator &&) = default; //!< Defaulted.
    ~trace_iterator() = default; //!< Defaulted.

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \param[in] matrix_iter The underlying matrix iterator.
     */
    constexpr trace_iterator(matrix_iter_t const matrix_iter) noexcept : matrix_iter{matrix_iter}
    {
        set_trace_direction(*matrix_iter);
    }

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \tparam other_matrix_iter_t The underlying matrix iterator type of `other`; the condition
     *                             `std::Constructible<matrix_iter_t, other_matrix_iter_t>` must evaluate to `true`.
     * \param[in] other The underlying matrix iterator.
     */
    template <TwoDimensionalMatrixIterator other_matrix_iter_t>
    //!\cond
        requires std::Constructible<matrix_iter_t, other_matrix_iter_t>
    //!\endcond
    constexpr trace_iterator(trace_iterator<other_matrix_iter_t> const other) noexcept :
        trace_iterator{other.matrix_iter}
    {}
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the current trace direction.
    reference operator*() const noexcept
    {
        return current_direction;
    }

    //!\brief Returns the current coordinate in two-dimensional space.
    constexpr matrix_coordinate coordinate() const noexcept
    {
        return matrix_iter.coordinate();
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advances the iterator by one.
    constexpr trace_iterator & operator++() noexcept
    {
        trace_directions old_dir = *matrix_iter;

        assert(old_dir != trace_directions::none);

        if (current_direction == trace_directions::up)
        {
            matrix_iter -= matrix_offset{row_index_type{1}, column_index_type{0}};
            // Set new trace direction if last position was up_open.
            if (static_cast<bool>(old_dir & trace_directions::up_open))
                set_trace_direction(*matrix_iter);
        }
        else if (current_direction == trace_directions::left)
        {
            matrix_iter -= matrix_offset{row_index_type{0}, column_index_type{1}};
            // Set new trace direction if last position was left_open.
            if (static_cast<bool>(old_dir & trace_directions::left_open))
                set_trace_direction(*matrix_iter);
        }
        else
        {
            assert(current_direction == trace_directions::diagonal);

            matrix_iter -= matrix_offset{row_index_type{1}, column_index_type{1}};
            set_trace_direction(*matrix_iter);
        }
        return *this;
    }

    //!\brief Returns an iterator advanced by one.
    constexpr trace_iterator operator++(int) noexcept
    {
        trace_iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Returns `true` if both iterators are equal, `false` otherwise.
    template <TwoDimensionalMatrixIterator other_matrix_iter_t>
    //!\cond
        requires std::Constructible<matrix_iter_t, other_matrix_iter_t> ||
                 std::Constructible<other_matrix_iter_t, matrix_iter_t>
    //!\endcond
    constexpr bool operator==(trace_iterator<other_matrix_iter_t> const & rhs) const noexcept
    {
        return *matrix_iter == *rhs.matrix_iter;
    }

    //!\brief Returns `true` if the pointed-to-element is seqan3::detail::trace_directions::none.
    constexpr bool operator==(std::ranges::default_sentinel_t const &) const noexcept
    {
        return *matrix_iter == trace_directions::none;
    }

    //!\brief copydoc operator==()
    constexpr friend bool operator==(std::ranges::default_sentinel_t const &, trace_iterator const & rhs) noexcept
    {
        return rhs == std::ranges::default_sentinel;
    }

    //!\brief Returns `true` if both iterators are not equal, `false` otherwise.
    template <TwoDimensionalMatrixIterator other_matrix_iter_t>
    //!\cond
        requires std::Constructible<matrix_iter_t, other_matrix_iter_t> ||
                 std::Constructible<other_matrix_iter_t, matrix_iter_t>
    //!\endcond
    constexpr bool operator!=(trace_iterator<other_matrix_iter_t> const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    //!\brief Returns `true` if the pointed-to-element is not seqan3::detail::trace_directions::none.
    constexpr bool operator!=(std::ranges::default_sentinel_t const &) const noexcept
    {
        return !(*this == std::ranges::default_sentinel);
    }

    //!\brief copydoc operator!=()
    constexpr friend bool operator!=(std::ranges::default_sentinel_t const &, trace_iterator const & rhs) noexcept
    {
        return !(rhs == std::ranges::default_sentinel);
    }
    //!\}

private:

    //!\brief Updates the current trace direction.
    void set_trace_direction(trace_directions const dir) noexcept
    {
        if (static_cast<bool>(dir & trace_directions::diagonal))
            current_direction = trace_directions::diagonal;
        else if (static_cast<bool>(dir & trace_directions::up) ||
                 static_cast<bool>(dir & trace_directions::up_open))
            current_direction = trace_directions::up;
        else if (static_cast<bool>(dir & trace_directions::left) ||
                 static_cast<bool>(dir & trace_directions::left_open))
            current_direction = trace_directions::left;
        else
            current_direction = trace_directions::none;
    }

    matrix_iter_t matrix_iter{}; //!< The underlying matrix iterator.
    trace_directions current_direction{};  //!< The current trace direction.
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Deduces the template argument from the passed iterator.
template <TwoDimensionalMatrixIterator matrix_iter_t>
trace_iterator(matrix_iter_t const) -> trace_iterator<matrix_iter_t>;
//!\}

} // namespace seqan3::detail
