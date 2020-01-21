// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::two_dimensional_matrix_iterator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>

namespace seqan3::detail
{

/*!\interface seqan3::detail::two_dimensional_matrix_iterator <>
 * \extends   std::random_access_iterator
 * \brief     A concept for iterators over a two dimensional matrix, e.g. seqan3::detail::two_dimensional_matrix
 * \ingroup   alignment_matrix
 *
 * This concept describes the requirements an iterator must fulfil in order to be used inside various parts
 * of the alignment algorithm, e.g. to compute the traceback path after filling the alignment matrix.
 */
/*!\name Requirements for seqan3::detail::two_dimensional_matrix_iterator
 * \brief You can expect these functions on all types that model seqan3::detail::two_dimensional_matrix_iterator.
 * \relates seqan3::detail::two_dimensional_matrix_iterator
 * \{
 */
/*!\fn iterator & operator+=(seqan3::detail::matrix_offset offset) noexcept
 * \brief Advances the iterator by `offset` in the respective dimension.
 *
 * \param[in] offset The matrix offset in a particular dimension.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn constexpr iterator operator+(seqan3::detail::matrix_offset offset) const noexcept
 * \brief Returns an iterator advanced by `offset` in the respective dimension.
 *
 * \param[in] offset The matrix offset in a particular dimension.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn constexpr friend iterator operator+(seqan3::detail::matrix_offset offset, iterator const iter) noexcept
 * \brief Returns an iterator advanced by `offset` in the respective dimension.
 *
 * \tparam iterator The two dimensional iterator type.
 * \param[in] offset The matrix offset in a particular dimension.
 * \param[in] iter   The iterator to advance by the offset.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn constexpr iterator & operator-=(seqan3::detail::matrix_offset offset) noexcept
 * \brief Advances the iterator by `offset` in the respective dimension.
 *
 * \param[in] offset The matrix offset in a particular dimension.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn constexpr iterator operator-(seqan3::detail::matrix_offset offset) const noexcept
 * \brief Returns an iterator advanced by `offset` in the respective dimension.
 *
 * \param[in] offset The matrix offset in a particular dimension.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn constexpr seqan3::detail::matrix_coordinate coordinate() const noexcept
 * \brief Returns the current position of the iterator as a two-dimensional matrix coordinate.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
//!\}
//!\cond
template <typename iter_t>
SEQAN3_CONCEPT two_dimensional_matrix_iterator =
    std::random_access_iterator<iter_t> &&
    requires(std::remove_reference_t<iter_t> it, std::remove_reference_t<iter_t> const cit, matrix_offset offset)
    {
        { it += offset };
        { it + offset };
        { offset + it };
        { cit + offset };
        { offset + cit };
        { it -= offset };
        { it - offset };
        { cit - offset };
        { it.coordinate() };
        { cit.coordinate() };

        requires std::same_as<decltype(it += offset), std::remove_reference_t<iter_t> &>;
        requires std::same_as<decltype(it + offset), std::remove_reference_t<iter_t>>;
        requires std::same_as<decltype(offset + it), std::remove_reference_t<iter_t>>;
        requires std::same_as<decltype(it -= offset), std::remove_reference_t<iter_t> &>;
        requires std::same_as<decltype(it - offset), std::remove_reference_t<iter_t>>;
        requires std::same_as<decltype(cit - offset), std::remove_reference_t<iter_t>>;
        requires std::same_as<decltype(it.coordinate()), matrix_coordinate>;
        requires std::same_as<decltype(cit.coordinate()), matrix_coordinate>;
    };
//!\endcond
} // namespace seqan3::detail
