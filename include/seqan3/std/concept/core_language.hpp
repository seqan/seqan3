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
 * \brief Adaptions of core language concepts from the STL / range-v3.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_concepts.hpp>
#include <range/v3/utility/functional.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup concept
 * \{
 */

/*!\brief Resolves to `ranges::Same<types...>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Same
 */
template <typename ...ts>
concept bool same_concept =                         static_cast<bool>(ranges::Same<ts...>());

/*!\brief Resolves to `ranges::DerivedFrom<type1, type2>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/DerivedFrom
 */
template <typename t, typename u>
concept bool derived_from_conept =                  static_cast<bool>(ranges::DerivedFrom<t, u>());

//!\brief Resolves to `ranges::ImplicitlyConvertibleTo<type1, type2>()`
template <typename t, typename u>
concept bool implicitly_convertible_to_concept =    static_cast<bool>(ranges::ImplicitlyConvertibleTo<t, u>());

//!\brief Resolves to `ranges::ExplicitlyConvertibleTo<type1, type2>()`
template <typename t, typename u>
concept bool explicitly_convertible_to_concept =    static_cast<bool>(ranges::ExplicitlyConvertibleTo<t, u>());

/*!\brief Subsumes, both, seqan3::implicitly_convertible_to_concept **and** seqan3::explicitly_convertible_to_concept.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/ConvertibleTo
 */
template <typename t, typename u>
concept bool convertible_to_concept =               implicitly_convertible_to_concept<t, u> &&
                                                    explicitly_convertible_to_concept<t, u>;

/*!\brief Resolves to `ranges::CommonReference<type1, type2, rest...>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/CommonReference
 */
template <typename t, typename u, typename ...rest>
concept bool common_reference_concept =             static_cast<bool>(ranges::CommonReference<t, u, rest...>());

/*!\brief Resolves to `ranges::Common<type1, type2, rest...>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Common
 */
template <typename t, typename u, typename ...rest>
concept bool common_concept =                       static_cast<bool>(ranges::Common<t, u, rest...>());

/*!\brief Resolves to `ranges::Integral<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Integral
 */
template <typename t>
concept bool integral_concept =                     static_cast<bool>(ranges::Integral<t>());

/*!\brief Resolves to `ranges::SignedIntegral<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/SignedIntegral
 */
template <typename t>
concept bool signed_integral_concept =              integral_concept<t> &&
                                                    static_cast<bool>(ranges::SignedIntegral<t>());

/*!\brief Resolves to `ranges::UnsignedIntegral<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/UnsignedIntegral
 */
template <typename t>
concept bool unsigned_integral_concept =            integral_concept<t> &&
                                                    static_cast<bool>(ranges::UnsignedIntegral<t>());

/*!\interface   seqan3::assignable_concept <>
 * \brief       Resolves to `std::is_assignable_v<t1, t2>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_assignable
 */
/*!\fn          t & operator=(t1 const & rhs)
 * \brief       Assignment operator.
 * \memberof    seqan3::assignable_concept
 * \param rhs   Right hand side parameter to assign.
 * \returns     Reference to self.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\interface   seqan3::trivially_assignable_concept <>
 * \extends     seqan3::assignable_concept
 * \brief       Resolves to `std::is_trivially_assignable_v<t1, t2>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_assignable
 */
/*!\interface   seqan3::nothrow_assignable_concept <>
 * \extends     seqan3::assignable_concept
 * \brief       Resolves to `std::is_nothrow_assignable_v<t1, t2>`.
 * \sa          http://en.cppreference.com/w/cpp/types/is_assignable
 */
//!\cond
template <typename t1, typename t2>
concept bool assignable_concept =                   std::is_assignable_v<t1, t2>;
template <typename t1, typename t2>
concept bool trivially_assignable_concept =         assignable_concept<t1, t2> &&
                                                    std::is_trivially_assignable_v<t1, t2>;
template <typename t1, typename t2>
concept bool nothrow_assignable_concept =           assignable_concept<t1, t2> &&
                                                    std::is_nothrow_assignable_v<t1, t2>;
//!\endcond

/*!\interface   seqan3::swappable_concept <>
 * \brief       Whether std::swap (or swap function in the associated namespace) works on a type.
 */
/*!\name Requirements for seqan3::swappable_concept
 * \brief You can expect these functions on all types that implement seqan3::swappable_concept.
 * \{
 */
/*!\fn          void std::swap(swappable_concept & lhs, swappable_concept & rhs)
 * \brief       Swaps the contents of two objects.
 * \relates     seqan3::swappable_concept
 * \param lhs   Left hand side parameter to swap.
 * \param rhs   Right hand side parameter to swap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
//!\cond
template <typename t>
concept bool swappable_concept =               requires (t && v1, t && v2) { std::swap(v1, v2); } ||
                                               requires (t && v1, t && v2) {      swap(v1, v2); };
//!\endcond

/*!\brief Resolves to `ranges::Swappable<type1, type2>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Swappable
 */
template <typename t1, typename t2>
concept bool swappable_with_concept =          swappable_concept<t1> &&
                                               swappable_concept<t2> &&
                                               common_reference_concept<std::remove_reference_t<t1> const &,
                                                                        std::remove_reference_t<t2> const &> &&
                                               requires (t1 && v1, t2 && v2)
                                               {
                                                   std::swap(v1, v2);
                                                   std::swap(v2, v1);

                                               } ||
                                               requires (t1 && v1, t2 && v2)
                                               {
                                                   swap(v1, v2);
                                                   swap(v2, v1);
                                               };
//!\}
}  // namespace seqan3
