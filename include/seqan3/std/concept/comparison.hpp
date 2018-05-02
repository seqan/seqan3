// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief  Concepts regarding a type's comparability with itself and other types.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concept/core_language.hpp>

namespace seqan3
{

/*!\addtogroup concept
 * \{
 */

/*!\interface   seqan3::weakly_equality_comparable_with_concept <>
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \brief       Requires the two operands to be comparable with `==` and `!=` in both directions.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/WeaklyEqualityComparableWith
 */
//!\cond
template <typename t1, typename t2>
concept bool weakly_equality_comparable_with_concept = requires (std::remove_reference_t<t1> const & v1,
                                                                 std::remove_reference_t<t2> const & v2)
{
    { v1 == v2 } -> bool &&;
    { v1 != v2 } -> bool &&;
    { v2 == v1 } -> bool &&;
    { v2 != v1 } -> bool &&;
};
//!\endcond

/*!\interface   seqan3::equality_comparable_concept <>
 * \brief       The same as seqan3::weakly_equality_comparable_with_concept<t,t>.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/EqualityComparable
 */
/*!\name Requirements for seqan3::equality_comparable_concept
 * \brief You can expect these functions on all types that implement seqan3::equality_comparable_concept.
 * \{
 */
/*!\fn          bool operator==(type const & lhs, type const & rhs)
 * \brief       (In-)Equality comparison.
 * \relates     seqan3::equality_comparable_concept
 * \param[in]   lhs Left hand side parameter to compare.
 * \param[in]   rhs Right hand side parameter to compare.
 * \returns     `true` or `false`, depending of the outcome of the comparison.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn      bool operator!=(type const & lhs, type const & rhs)
 * \relates seqan3::equality_comparable_concept
 * \copydoc seqan3::equality_comparable_concept::operator==
 */
//!\}
//!\cond
template <typename t>
concept bool equality_comparable_concept = weakly_equality_comparable_with_concept<t, t>;
//!\endcond

/*!\interface   seqan3::equality_comparable_with_concept <>
 * \extends     seqan3::weakly_equality_comparable_with_concept
 * \brief       Requires seqan3::weakly_equality_comparable_with_concept<t1,t2>, but also that t1 and t2, as well as
 *              their common_reference_t satisfy seqan3::equality_comparable_concept.
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/EqualityComparable
 */
//!\cond
template <typename t1, typename t2>
concept bool equality_comparable_with_concept =     equality_comparable_concept<t1> &&
                                                    equality_comparable_concept<t2> &&
                                                    common_reference_concept<std::remove_reference_t<t1> const &,
                                                                             std::remove_reference_t<t2> const &> &&
                                                    equality_comparable_concept<
                                                        ranges::common_reference_t<std::remove_reference_t<t1> const &,
                                                                                   std::remove_reference_t<t2> const &>> &&
                                                    weakly_equality_comparable_with_concept<t1, t2>;
//!\endcond

} // namespace seqan3

namespace seqan3::detail
{

//!\cond
template <typename t1, typename t2>
concept bool weakly_ordered_with_concept =          requires (std::remove_reference_t<t1> const & v1,
                                                              std::remove_reference_t<t2> const & v2)
{
    { v1 <  v2 } -> bool &&;
    { v1 <= v2 } -> bool &&;
    { v2 >  v1 } -> bool &&;
    { v2 >= v1 } -> bool &&;
};
//!\endcond

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface   seqan3::strict_totally_ordered_concept
 * \extends     seqan3::equality_comparable_concept
 * \brief       Requires weak order and seqan3::equality_comparable_concept.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/StrictTotallyOrdered
 */
/*!\name Requirements for seqan3::strict_totally_ordered_concept
 * \brief You can expect these functions on all types that implement seqan3::strict_totally_ordered_concept.
 * \{
 */
/*!\fn          bool operator<(type const & lhs, type const & rhs)
 * \brief       Less-than, greater-than and -or-equal comparisons.
 * \relates     seqan3::strict_totally_ordered_concept
 * \param[in]   lhs Left hand side parameter to compare.
 * \param[in]   rhs Right hand side parameter to compare.
 * \returns     `true` or `false`, depending of the outcome of the comparison.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn      bool operator<=(type const & lhs, type const & rhs)
 * \relates seqan3::strict_totally_ordered_concept
 * \copydoc seqan3::strict_totally_ordered_concept::operator<
 */
/*!\fn      bool operator>(type const & lhs, type const & rhs)
 * \relates seqan3::strict_totally_ordered_concept
 * \copydoc seqan3::strict_totally_ordered_concept::operator<
 */
/*!\fn      bool operator>=(type const & lhs, type const & rhs)
 * \relates seqan3::strict_totally_ordered_concept
 * \copydoc seqan3::strict_totally_ordered_concept::operator<
 */
//!\}
//!\cond
template <typename t>
concept bool strict_totally_ordered_concept =        equality_comparable_concept<t> &&
                                                     detail::weakly_ordered_with_concept<t, t>;
//!\endcond

/*!\interface   seqan3::totally_ordered_with_concept
 * \extends     seqan3::equality_comparable_with_concept
 * \brief       Requires seqan3::weakly_ordered_concept and seqan3::equality_comparable_concept.
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/StrictTotallyOrdered
 */
//!\cond
template <typename t1, typename t2>
concept bool strict_totally_ordered_with_concept =   strict_totally_ordered_concept<t1> &&
                                                     strict_totally_ordered_concept<t2> &&
                                                     strict_totally_ordered_concept<
                                                        ranges::common_reference_t<std::remove_reference_t<t1> const &,
                                                                                   std::remove_reference_t<t2> const &>> &&
                                                     detail::weakly_ordered_with_concept<t1, t2> &&
                                                     equality_comparable_with_concept<t1, t2>;
//!\endcond

//!\}

}  // namespace seqan3
