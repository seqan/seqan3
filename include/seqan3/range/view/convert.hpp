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
 * \brief Provides seqan3::view::convert.
 */

#pragma once

#include <range/v3/view/transform.hpp>

#include <seqan3/core/concept/core.hpp>

namespace seqan3::view
{

/*!\brief A view that converts each element in the input range (implicitly or via `static_cast`).
 * \tparam irng_t The type of the range being processed. See below for requirements. [template parameter is omitted in pipe notation]
 * \param irange The range being processed. [parameter is omitted in pipe notation]
 * \returns A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * \par View properties
 *
 * |                     | `irng_t` (range input type)   | `rrng_t` (range return type)                              |
 * |---------------------|-------------------------------|-----------------------------------------------------------|
 * | range               | seqan3::input_range_concept   | seqan3::view_concept + all range concepts met by `irng_t` |
 * | `range_reference_t` | *convertible to* `out_t`      | `out_t`                                                   |
 *
 * * The input properties are **requirements** on the range input type.
 * * The return properties are **guarantees** given on the range return type.
 * * for more details, see \ref view.
 *
 * \par Example
 *
 * Convert from `int` to `bool`:
 * ```cpp
 *   // convert from int to bool
 *   std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};
 *
 *   // pipe notation
 *   auto v = vec | view::convert<bool>; // == [1, 1, 0, 1, 0, 0, 1, 1, 1];
 *
 *   // function notation and immediate conversion to vector again
 *   std::vector<bool> v2(view::convert<bool>(vec));
 *
 *   // combinability
 *   auto v3 = vec | view::convert<bool> | ranges::view::reverse; // == [1, 1, 1, 0, 0, 1, 0, 1, 1];
 * ```
 *
 * Convert from seqan3::nucl16 to seqan3::dna5:
 * ```cpp
 *   nucl16_vector vec2{"ACYGTN"_nucl16};
 *   auto v4 = vec2 | view::convert<dna5>; // == "ACNGTN"_dna5
 * ```
 * \hideinitializer
 */
template <typename out_t>
auto const convert = ranges::view::transform([] (auto const & in) -> out_t
{
    if constexpr (implicitly_convertible_to_concept<std::remove_reference_t<decltype(in)>, out_t>)
        return in;
    else
        return static_cast<out_t>(in);
});

} // namespace seqan3::view
