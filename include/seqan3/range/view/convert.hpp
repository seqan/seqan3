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
 * \brief Provides seqan3::view::convert.
 */

#pragma once

<<<<<<< HEAD
#include <range/v3/view/transform.hpp>

#include <seqan3/std/concept/core_language.hpp>
=======
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/view/transform.hpp>
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view that converts each element in the input range (implicitly or via `static_cast`).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * ### View properties
 *
<<<<<<< HEAD
 * | range concepts and reference_t      | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |-------------------------------------|:-------------------------------------:|:-------------------------------:|
 * | seqan3::input_range_concept         | *required*                            | *preserved*                     |
 * | seqan3::forward_range_concept       |                                       | *preserved*                     |
 * | seqan3::bidirectional_range_concept |                                       | *preserved*                     |
 * | seqan3::random_access_range_concept |                                       | *preserved*                     |
 * |                                     |                                       |                                 |
 * | seqan3::view_concept                |                                       | *guaranteed*                    |
 * | seqan3::sized_range_concept         |                                       | *preserved*                     |
 * | seqan3::bounded_range_concept       |                                       | *preserved*                     |
 * | seqan3::output_range_concept        |                                       | *lost*                          |
 * | seqan3::const_iterable_concept      |                                       | *preserved*                     |
 * |                                     |                                       |                                 |
 * | seqan3::reference_t                 | seqan3::convertible_to_concept<out_t> | `out_t`                         |
=======
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |---------------------------------|:-------------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                     |
 * | std::ranges::ForwardRange       |                                       | *preserved*                     |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                     |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | std::ranges::View               |                                       | *guaranteed*                    |
 * | std::ranges::SizedRange         |                                       | *preserved*                     |
 * | std::ranges::CommonRange        |                                       | *preserved*                     |
 * | std::ranges::OutputRange        |                                       | *lost*                          |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | seqan3::reference_t             | seqan3::convertible_to_concept<out_t> | `out_t`                         |
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * Convert from `int` to `bool`:
 * \snippet test/snippet/range/view/convert.cpp int_to_bool
 *
 * Convert from seqan3::dna15 to seqan3::dna5:
 * \snippet test/snippet/range/view/convert.cpp 15_to_5
 * \hideinitializer
 */
template <typename out_t>
auto const convert = view::transform([] (auto const & in) -> out_t
{
    if constexpr (implicitly_convertible_to_concept<std::remove_reference_t<decltype(in)>, out_t>)
        return in;
    else
        return static_cast<out_t>(in);
});

//!\}

} // namespace seqan3::view
