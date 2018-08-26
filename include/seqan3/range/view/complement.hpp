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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::complement.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that converts a range of nucleotides to their complement.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * Calls seqan3::nucleotide_concept::complement() on every element of the input range.
 *
 * ### View properties
 *
 * This view is a **deep view:** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *preserved*                                        |
 * | std::ranges::OutputRange        |                                       | *lost*                                             |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             | seqan3::nucleotide_concept            | std::remove_reference_t<seqan3::reference_t<urng_t>> |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * ```cpp
 *  dna5_vector foo{"ACGTA"_dna5};
 *
 *  // pipe notation
 *  auto v = foo | view::complement;                                  // == "TGCAT"
 *
 *  // function notation
 *  dna5_vector v2(view::complement(foo));                            // == "TGCAT"
 *
 *  // generate the reverse complement:
 *  dna5_vector v3 = foo | view::complement | view::reverse;  // == "TACGT"
 * ```
 * \hideinitializer
 */

inline auto const complement = deep{view::transform([] (auto && in)
{
    static_assert(nucleotide_concept<std::remove_const_t<decltype(in)>>,
                  "The innermost value type must satisfy the nucleotide_concept.");
    // call element-wise complement from the nucleotide_concept
    using seqan3::complement;
    return complement(in);
})};

//!\}

} // namespace seqan3::view
