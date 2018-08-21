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

#include <seqan3/std/concept/comparison.hpp>

namespace seqan3
{

/*!\addtogroup concept
 * \{
 */

/*!\interface   seqan3::invocable_concept <>
 * \brief       Specifies whether the given callable is invocable with the given arguments.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Invocable
 */
//!\cond
template <typename fun_t, typename ...arg_types>
concept bool invocable_concept = requires(fun_t && f, arg_types &&... args)
{
    ranges::invoke(std::forward<fun_t>(f), std::forward<arg_types>(args)...);
      /* not required to be equality preserving */
};
//!\endcond

/*!\interface   seqan3::regular_invocable_concept <>
 * \extends     seqan3::invocable_concept
 * \brief       Specifies whether the given callable is invocable with the given arguments and equality preserving
 *              (invocations change neither the callable, nor the arguments).
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/RegularInvocable
 */
//!\cond
template <typename fun_t, typename ...arg_types>
concept bool regular_invocable_concept = invocable_concept<fun_t, arg_types...> &&
                                         requires(fun_t && f, arg_types &&... args)
{
    ranges::invoke(static_cast<std::remove_reference_t<fun_t> const &>(f), std::forward<arg_types>(args)...);
      /* required to be equality preserving */
};
//!\endcond

/*!\interface   seqan3::predicate_concept <>
 * \extends     seqan3::regular_invocable_concept
 * \brief       Specifies whether the given callable is seqan3::regular_invocable_concept and returns bool.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Predicate
 */
//!\cond
template <typename fun_t, typename ...arg_types>
concept bool predicate_concept = regular_invocable_concept<fun_t, arg_types...> &&
                                 boolean_concept<std::result_of_t<fun_t&&(arg_types&&...)>>;
//!\endcond

/*!\interface   seqan3::relation_concept <>
 * \extends     seqan3::predicate_concept
 * \brief       Specifies that fun_t defines a binary relation over the set of expressions whose type and value category
 *              are those encoded by either t or u.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Relation
 */
//!\cond
template <typename fun_t, typename t, typename u>
concept bool relation_concept = predicate_concept<fun_t, t, t> &&
                                predicate_concept<fun_t, u, u> &&
                                common_reference_concept<std::remove_reference_t<t> const &,
                                                         std::remove_reference_t<u> const &> &&
                                predicate_concept<fun_t,
                                    ranges::common_reference_t<std::remove_reference_t<t> const &,
                                                               std::remove_reference_t<u> const &>,
                                    ranges::common_reference_t<std::remove_reference_t<t> const &,
                                                               std::remove_reference_t<u> const &>> &&
                                predicate_concept<fun_t, t, u> &&
                                predicate_concept<fun_t, u, t>;
//!\endcond
//!\}
}  // namespace seqan3
