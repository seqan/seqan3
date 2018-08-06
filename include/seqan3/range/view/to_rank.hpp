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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::to_rank.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that calls seqan3::to_rank() on each element in the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * ### View properties
 *
 * This view is a **deep view:** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | range concepts and reference_t      | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |-------------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | seqan3::input_range_concept         | *required*                            | *preserved*                                        |
 * | seqan3::forward_range_concept       |                                       | *preserved*                                        |
 * | seqan3::bidirectional_range_concept |                                       | *preserved*                                        |
 * | seqan3::random_access_range_concept |                                       | *preserved*                                        |
 * | seqan3::contiguous_range_concept    |                                       | *lost*                                             |
 * |                                     |                                       |                                                    |
 * | seqan3::viewable_range_concept      | *required*                            | *guaranteed*                                       |
 * | seqan3::view_concept                |                                       | *guaranteed*                                       |
 * | seqan3::sized_range_concept         |                                       | *preserved*                                        |
 * | seqan3::common_range_concept        |                                       | *preserved*                                        |
 * | seqan3::output_range_concept        |                                       | *lost*                                             |
 * | seqan3::const_iterable_concept      |                                       | *preserved*                                        |
 * |                                     |                                       |                                                    |
 * | seqan3::reference_t                 | seqan3::alphabet_concept              | seqan3::underlying_rank_t<seqan3::value_type_t<urng_t>> |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * \par Example
 * ```cpp
 * dna4_vector vec = "ACTTTGATA"_dna4;
 * auto v = vec | view::to_rank | view::convert<unsigned>;
 * std::cout << v << '\n'; // [0,1,3,3,3,2,0,3,0]
 *
 * std::vector<illumina18> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
 * auto v3 = qvec | view::to_rank | view::convert<unsigned>;
 * std::cout << v3 << '\n'; // [0,7,5,3,7,4,30,16,23]
 *
 * std::vector<dna4q> qcvec{{dna4::C, 0}, {dna4::A, 7}, {dna4::G, 5}, {dna4::T, 3}, {dna4::G, 7}, {dna4::A, 4}, {dna4::C, 30}, {dna4::T, 16}, {dna4::A, 23}};
 * auto v4 = qcvec | view::to_rank | view::convert<unsigned>;
 * std::cout << v4 << '\n'; // [1,28,22,15,30,16,121,67,92]
 * ```
 * We also convert to unsigned here, because the seqan3::underlying_rank_t is often `uint8_t` which is
 * often implemented as `unsigned char` and thus will not be printed as a number by default.
 * \hideinitializer
 */
inline auto const to_rank = deep{view::transform([] (alphabet_concept const in)
{
    return seqan3::to_rank(in);
})};

//!\}

} // namespace seqan3::view
