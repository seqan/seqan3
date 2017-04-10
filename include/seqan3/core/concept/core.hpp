// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

/*!\file core/concept/core.hpp
 * \brief Adaptions of core concepts from the Ranges TS.
 * \ingroup core
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_concepts.hpp>
#include <range/v3/utility/functional.hpp>

namespace seqan3
{

//!\name Core Language Concepts
//!\{
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

/*!\brief Resolves to `ranges::ConvertibleTo<type1, type2>()` 
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

/*!\brief Resolves to `ranges::Assignable<type1, type2>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Assignable
 */
template <typename t, typename u>
concept bool assignable_concept =                   static_cast<bool>(ranges::Assignable<t, u>());

/*!\brief Resolves to `ranges::Swappable<type1, type2>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Swappable
 */
template <typename t, typename u = t>
concept bool swappable_concept =                    static_cast<bool>(ranges::Swappable<t, u>());

//!\}  // Core Language Concepts.

//!\name Comparison concepts.
//!\{

/*!\brief Resolves to `ranges::WeaklyEqualityComparable<type1, type2>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/WeaklyEqualityComparable
 */
template <typename t, typename u>
concept bool weakly_equality_comparable_concept =   static_cast<bool>(ranges::WeaklyEqualityComparable<t, u>());

/*!\brief Resolves to `ranges::EqualityComparable<type1, type2>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/EqualityComparable
 */
template <typename t, typename u = t>
concept bool equality_comparable_concept =          static_cast<bool>(ranges::EqualityComparable<t, u>());

/*!\brief Resolves to `ranges::WeaklyOrdered<type1, type2>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/WeaklyOrdered
 */
template <typename t, typename u = t>
concept bool weakly_ordered_concept =               static_cast<bool>(ranges::WeaklyOrdered<t>());

/*!\brief Resolves to `ranges::TotallyOrdered<type1, type2>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/TotallyOrdered
 */
template <typename t, typename u = t>
concept bool totally_ordered_concept =              equality_comparable_concept<t, u> &&
                                                    weakly_ordered_concept<t, u> &&
                                                    static_cast<bool>(ranges::TotallyOrdered<t, u>());
//!\}  // Comparison Concepts.

//!\name Object Concepts.
//!\{

/*!\brief Resolves to `ranges::Destructible<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Destructible
 */
template <typename t>
concept bool destructible_concept =                 static_cast<bool>(ranges::Destructible<t>());

/*!\brief Resolves to `ranges::Constructible<type, args...>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Constructible
 */
template <typename t, typename ...args>
concept bool constructible_concept =                destructible_concept<t> &&
                                                    static_cast<bool>(ranges::Constructible<t, args...>());

/*!\brief Resolves to `ranges::DefaultConstructible<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/DefaultConstructible
 */
template <typename t>
concept bool default_constructible_concept =        constructible_concept<t> &&
                                                    static_cast<bool>(ranges::DefaultConstructible<t>());

/*!\brief Resolves to `ranges::MoveConstructible<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/MoveConstructible
 */
template <typename t>
concept bool move_constructible_concept =           static_cast<bool>(ranges::MoveConstructible<t>());

/*!\brief Resolves to `ranges::CopyConstructible<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/CopyConstructible
 */
template <typename t>
concept bool copy_constructible_concept =           move_constructible_concept<t> &&
                                                    static_cast<bool>(ranges::CopyConstructible<t>());

/*!\brief Resolves to `ranges::Movable<types...>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Movable
 */
template <typename t>
concept bool movable_concept =                      move_constructible_concept<t> &&
                                                    static_cast<bool>(ranges::Movable<t>());

/*!\brief Resolves to `ranges::Copyable<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Copyable
 */
template <typename t>
concept bool copyable_concept =                     movable_concept<t> &&
                                                    copy_constructible_concept<t> &&
                                                    static_cast<bool>(ranges::Copyable<t>());

/*!\brief Resolves to `ranges::SemiRegular<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/SemiRegular
 */
template <typename t>
concept bool semi_regular_concept =                 copyable_concept<t> &&
                                                    default_constructible_concept<t>;

/*!\brief Resolves to `ranges::Regular<type>()` 
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Regular
 */
template <typename t>
concept bool regular_concept =                      semi_regular_concept<t> &&
                                                    equality_comparable_concept<t>;
//!\}  Object Concepts.

//!\name Callable Concepts.
//!\{

/*!\brief Resolves to `ranges::Invocable<func, ...args>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Invocable
 */
template <typename f, typename ...args>
concept bool invocable_concept =                    static_cast<bool>(ranges::Invocable<f, args...>());

/*!\brief Resolves to `ranges::RegularInvocable<func, ...args>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/RegularInvocable
 */
template <typename f, typename ...args>
concept bool regular_invocable_concept =            invocable_concept<f, args...> &&
                                                    static_cast<bool>(ranges::RegularInvocable<f, args...>());

/*!\brief Resolves to `ranges::Predicate<func, ...args>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Predicate
 */
template <typename f, typename ...args>
concept bool predicate_concept =                    regular_invocable_concept<f, args...> &&
                                                    static_cast<bool>(ranges::Predicate<f, args...>());

/*!\brief Resolves to `ranges::Relation<func, type1, type2>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Relation
 */
template <typename f, typename t, typename u>
concept bool relation_concept =                     static_cast<bool>(ranges::Relation<f, t, u>());

//!\}  Callable Concepts.

}  // namespace seqan3

#ifndef NDEBUG

#include <seqan3/core/concepts/core_detail.hpp>

#endif // NDEBUG
