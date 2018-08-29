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
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::view::get.
 */

#pragma once

#include <seqan3/std/view/transform.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/concept/tuple.hpp>

namespace seqan3::view
{
/*!\name General purpose views
 * \{
 */

/*!\brief               A view calling std::get on each element in a range.
 * \tparam size_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam index        The index to get.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of elements where every element is the result of calling std::get<index> on the underlying element.
                        See below for the properties of the returned range.
 * \ingroup view
 *
 * ### View properties
 *
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                          |
 * |---------------------------------|:-------------------------------------:|:-------------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                             |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                             |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                             |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                             |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                                  |
 * |                                 |                                       |                                                         |
 * | std::ranges::ViewableRange      | *required*                            | *preserved*                                             |
 * | std::ranges::View               |                                       | *preserved*                                             |
 * | std::ranges::SizedRange         |                                       | *preserved*                                             |
 * | std::ranges::CommonRange        |                                       | *preserved*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                             |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                             |
 * |                                 |                                       |                                                         |
 * | seqan3::reference_t             | seqan3::tuple_like_concept            | std::tuple_element_t<index, seqan3::reference_t<urng_t>>|
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \snippet test/snippet/range/view/get.cpp usage
 * \hideinitializer
 */
template <size_t index>
inline auto const get = view::transform([] (auto && in)
    -> std::conditional_t<std::is_lvalue_reference_v<decltype(in)> && !std::is_const_v<decltype(in)>,
                          std::tuple_element_t<index, remove_cvref_t<decltype(in)>> &,
                          std::tuple_element_t<index, remove_cvref_t<decltype(in)>>>
{
    using std::get;
    static_assert(tuple_like_concept<decltype(in)>,
                  "You may only pass ranges to view::get whose reference_t models the tuple_like_concept.");
    return std::get<index>(in);
});

//!\}

} // namespace seqan3::view
