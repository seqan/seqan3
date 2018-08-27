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
 * \brief Provides seqan3::view::take_exactly and seqan3::view::take_exactly_or_throw.
 */

#pragma once

#include <seqan3/range/view/take.hpp>

// ============================================================================
//  view::take_exactly (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns the first `size` elements from the underlying range (or less if the
 *                      underlying range is shorter); also provides size information.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] size      The target size of the view.
 * \returns             Up to `size` elements of the underlying range.
 * \ingroup view
 *
 * \details
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | ***guaranteed***                                   |
 * | std::ranges::CommonRange        |                                       | *lost*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * The difference to seqan3::view::take is that this view always exposes size information (while
 * seqan3::view::take never does so). You should only use this if you know that the underlying range will always be
 * at least `size` long.
 *
 * For seqan3::view::take_exactly if the underlying range is shorter than `size`, the behaviour is undefined.
 * seqan3::view::take_exactly_or_throw is a safer alternative, because it throws an exception when an iterator before
 * the `size`-th one compares equal to the end sentinel; and it also throws on construction if it knows that the
 * underlying range is smaller.
 *
 * ### Example
 *
 * ```cpp
 * std::string vec{"foobar"};
 * auto v = vec | view::take_exactly(3);        // or view::take_exactly_or_throw
 * std::cout << v << '\n';                      // [f,o,o]
 * std::cout << ranges::size(v) << v << '\n';   // 3
 * ```
 *
 * The behaviour differs when the underlying sequence is shorter:
 *
 * ```cpp
 * std::string vec{"foo"};
 * auto v = vec | view::take_exactly(4);
 * std::cout << v << '\n';                          // [f,o,o]
 * std::cout << ranges::size(v) << v << '\n';       // 4 <- here be dragons!
 *
 * auto v2 = vec | view::take_exactly_or_throw(4);  // throws immediately on construction
 * ```
 *
 * \hideinitializer
 */
inline auto constexpr take_exactly = detail::take_fn<true, false>{};

//!\}

// ============================================================================
//  view::take_exactly_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns the first `size` elements from the underlying range and also exposes size
 *        information; throws if the underlying range is smaller than `size`.
 * \throws seqan3::unexpected_end_of_input If the underlying range is smaller than `size`.
 * \ingroup view
 *
 * \copydetails seqan3::view::take_exactly
 * \hideinitializer
 */
inline auto constexpr take_exactly_or_throw = detail::take_fn<true, true>{};

//!\}
} // namespace seqan3::view
