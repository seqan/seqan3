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

/*!\file
 * \ingroup view
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::rank_to.
 */

#pragma once

#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3::view
{

/*!\brief A view over an alphabet, given a range of ranks.
 * \tparam alphabet_type The type of the desired alphabet, must satisfy seqan3::alphabet_concept.
 * \param input_range The range you wish to convert, elements be must convertible to alphabet_type's
 * seqan3::underlying_rank_t.
 * \returns A view over alphabet_type, created from it's rank representation.
 * \details
 * \par View properties
 * * view type: same input_range
 * * value type: alphabet_type
 * * `const` iterable: yes
 * \par Complexity
 * Linear in the size if the input range (\f$O(n)\f$).
 * \par Exceptions
 * Strong exception guarantee (does not modify data).
 * \par Thread safety
 * Does not modify data.
 * \par Example
 * ```cpp
 * std::vector<int> vec{0, 1, 3, 3, 3, 2, 0, 3, 0};
 * auto v1 = vec | view::rank_to<dna4>; // == "ACTTTGATA"_dna4
 * auto v2 = vec | view::rank_to<dna5>; // == "ACTTTGATA"_dna5
 * ```
 * \hideinitializer
 */
template <typename alphabet_type>
//!\cond
    requires alphabet_concept<alphabet_type>
//!\endcond
auto const rank_to = ranges::view::transform([] (underlying_rank_t<alphabet_type> const in) -> alphabet_type
{
    return assign_rank(alphabet_type{}, in);
});

} // namespace seqan3::view
