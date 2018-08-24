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
 * \brief Provides seqan3::view::reverse.
 */

#pragma once

#include <range/v3/view/reverse.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::view
{

/*!\brief A range adaptor that presents the underlying range in reverse order.
 * \tparam           urng_t The type of the range being processed. See below for requirements. [template parameter
 *                          is omitted in pipe notation]
 * \param[in]        urange The range being processed. [parameter is omitted in pipe notation]
 * \returns A view of the elements of the underlying range in reverse order.
 * \ingroup core_view
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                       |
 * |---------------------------------|:-------------------------------------:|:----------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                          |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                                          |
 * | std::ranges::BidirectionalRange | *required*                            | *preserved*                                          |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                          |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                               |
 * |                                 |                                       |                                                      |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                         |
 * | std::ranges::View               |                                       | *guaranteed*                                         |
 * | std::ranges::SizedRange         |                                       | *preserved*                                          |
 * | std::ranges::CommonRange        |                                       | *preserved*                                          |
 * | std::ranges::OutputRange        |                                       | *lost*                                               |
 * | seqan3::const_iterable_concept  |                                       | *preserved if std::ranges::CommonRange<urng_t>*  |
 * |                                 |                                       |                                                      |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                          |
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
 * std::string s{"ACTNTGATAN"};
 * auto v1 = s | view::reverse; // == "NATAGTNTCA"
 * ```
 * \hideinitializer
 */

inline constexpr auto reverse = ranges::view::reverse;

} // namespace seqan3::view
