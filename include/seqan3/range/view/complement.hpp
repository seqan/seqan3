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
 * \brief Provides seqan3::view::complement.
 */

#pragma once

#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace seqan3::view
{

/*!\brief A view that converts a range of nucleotides to their complement.
 * \param input_range The range you wish to convert, must satisfy seqan3::input_range_concept and the value_type
 * must satisfy seqan3::nucleotide_concept.
 * \ingroup view
 * \details
 * \par View properties
 * * view type: same as input_range
 * * value type: remove_reference_t<value_type_t<input_range>>
 * * `const` iterable: yes
 * \par Complexity
 * Linear in the size if the input range (\f$O(n)\f$).
 * \par Exceptions
 * Strong exception guarantee (does not modify data).
 * \par Thread safety
 * Does not modify data.
 * \par Example
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
 *  dna5_vector v3 = foo | view::complement | ranges::view::reverse;  // == "TACGT"
 * ```
 * \hideinitializer
 */

//TODO enforce nucleotide_concept via c++2a and/or terse concept syntax
auto const complement = ranges::view::transform([] (auto const & in)
{
    using seqan3::complement;
    return complement(in);
});

} // namespace seqan3::view
