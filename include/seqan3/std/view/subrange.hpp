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
 * \brief Provides seqan3::view::subrange.
 */

#pragma once

#include <range/v3/iterator_range.hpp>

#include <seqan3/std/iterator>

namespace seqan3::view
{

/*!\brief Create a view from a pair of iterator and sentinel.
 * \tparam   it_t Type of the iterator; must satisfy std::Iterator.
 * \tparam  sen_t Type of the sentinel; must satisfy std::Sentinel with it_t.
 * \param[in]  it The iterator on the underlying range.
 * \param[in] sen The sentinel on the underlying range
 * \returns  A view of the elements between it_t and sen_t.
 * \ingroup core_view
 *
 * ### View properties
 *
 * This view is **source-only**, it can only be at the beginning of a pipe of range transformations.
 *
 * | range concepts and reference_t  | `rrng_t` (returned range type)                     |
 * |---------------------------------|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *preserved*                                        |
 * | std::ranges::ForwardRange       | *preserved*                                        |
 * | std::ranges::BidirectionalRange | *preserved*                                        |
 * | std::ranges::RandomAccessRange  | *preserved*                                        |
 * | std::ranges::ContiguousRange    | *preserved*                                        |
 * |                                 |                                                    |
 * | std::ranges::ViewableRange      | *guaranteed*                                       |
 * | std::ranges::View               | *guaranteed*                                       |
 * | std::ranges::SizedRange         | *preserved*                                        |
 * | std::ranges::CommonRange        | *preserved*                                        |
 * | std::ranges::OutputRange        | *preserved*                                        |
 * | seqan3::const_iterable_concept  | *preserved*                                        |
 * |                                 |                                                    |
 * | seqan3::reference_t             | seqan3::value_type_t<it_t>                         |
 *
 * Preservation in this table refers to the properties of the iterator/sentinel pair.
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### STD module
 *
 * This entity will likely be part of C++20 or C++23. It's API will track current proposals and not be stable
 * within SeqAn3 releases. It is implemented via the range-v3 library or the standard library (if available).
 * Should it become clear that it will not become part of a future standard, it will migrate to a regular SeqAn3
 * module.
 *
 * ### Example
 *
 * ```cpp
 * dna4_vector s{"ACTTTGATAN"_dna4};
 * auto v1 = subrange{begin(s) + 2, end(s)} | view::to_char; // == "TTTGATAA"
 * ```
 * \hideinitializer
 */
template <std::Iterator it_t, std::Sentinel<it_t> sen_t>
using subrange = ranges::iterator_range<it_t, sen_t>;
//TODO change to ranges::subrange once that has arrived

} // namespace seqan3::view
