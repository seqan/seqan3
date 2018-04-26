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
 * \brief Adaptions of core concepts from the Ranges TS.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_concepts.hpp>
#include <range/v3/utility/functional.hpp>

namespace seqan3
{

/*!\addtogroup concept
 * \{
 */

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

//!\}
}  // namespace seqan3
