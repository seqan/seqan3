// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/trace_iterator_base.hpp>

namespace seqan3::detail
{

/*!\brief A trace iterator an unbanded trace matrix.
 * \ingroup alignment_matrix
 *
 * \tparam matrix_iter_t The wrapped matrix iterator; must model seqan3::detail::two_dimensional_matrix_iterator and
 *                       the iterator's value type must be the same as seqan3::detail::trace_directions, i.e.
 *                       `std::same_as<std::iter_value_t<matrix_iter_t>, trace_directions>` must evaluate to `true`.
 *
 * \details
 *
 * This iterator follows a given trace in an unbanded trace matrix.
 */
template <two_dimensional_matrix_iterator matrix_iter_t>
class trace_iterator : public trace_iterator_base<trace_iterator<matrix_iter_t>, matrix_iter_t>
{
private:
    static_assert(std::same_as<std::iter_value_t<matrix_iter_t>, trace_directions>,
                  "Value type of the underlying iterator must be trace_directions.");

    //!\brief The type of the base class.
    using base_t = trace_iterator_base<trace_iterator<matrix_iter_t>, matrix_iter_t>;

public:
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
    explicit constexpr trace_iterator(matrix_iter_t const matrix_iter) noexcept : base_t{matrix_iter}
    {}

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \tparam other_matrix_iter_t The underlying matrix iterator type of `other`; the condition
     *                             `std::constructible_from<matrix_iter_t, other_matrix_iter_t>` must evaluate to `true`.
     * \param[in] other The underlying matrix iterator.
     *
     * \details
     *
     * Allows the conversion of non-const to const iterator.
     */
    template <two_dimensional_matrix_iterator other_matrix_iter_t>
    //!\cond
        requires std::constructible_from<matrix_iter_t, other_matrix_iter_t>
    //!\endcond
    constexpr trace_iterator(trace_iterator<other_matrix_iter_t> const other) noexcept : base_t{other}
    {}
    //!\}
};

} // namespace seqan3::detail
