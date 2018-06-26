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
 * \brief Provides seqan3::tuple_like_concept.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/std/concept/comparison.hpp>
#include <seqan3/std/concept/core_language.hpp>
#include <seqan3/std/concept/object.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// tuple_like_concept
// ----------------------------------------------------------------------------

/*!\interface   seqan3::tuple_like_concept
 * \extends     seqan3::strict_totally_ordered_concept
 * \ingroup     core
 * \brief       Whether a type behaves like a tuple.
 *
 * \details
 *
 * Types that meet this concept are std::tuple, std::pair, std::array, seqan3::pod_tuple, seqan3::record.
 */
/*!\name Requirements for seqan3::tuple_like_concept
 * \brief You can expect these (meta-)functions on all types that implement seqan3::tuple_like_concept.
 * \{
 */
/*!\var         size_t std::tuple_size_v<type>
 * \brief       A unary type trait that holds the number of elements in the tuple.
 * \tparam      type The tuple-like type.
 * \relates     seqan3::tuple_like_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\typedef     std::tuple_elment_t<i, type>
 * \brief       A transformation trait that holds the type of elements in the tuple.
 * \tparam      i Index of the queried element type.
 * \tparam      type The tuple-like type.
 * \relates     seqan3::tuple_like_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 * \attention This constraint is **not enforced** since empty tuples are valid.
 */
/*!\fn              auto && std::get<i>(type && val)
 * \brief           Return the i-th element of the tuple.
 * \relates         seqan3::tuple_like_concept
 * \tparam          i The index of the element to return (of type `size_t`).
 * \param[in,out]   val The tuple-like object to operate on.
 * \returns         The i-th value in the tuple.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 * \attention This constraint is **not enforced** since empty tuples are valid.
 */
//!\}
//!\cond
template <typename t>
concept bool tuple_like_concept = strict_totally_ordered_concept<t> &&
                                  requires (t v)
{
    { std::tuple_size_v<std::remove_reference_t<t>> } -> size_t;
};
//!\endcond

} // namespace seqan3
