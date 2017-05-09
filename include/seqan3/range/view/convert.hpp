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

/*!\file range/view/convert.hpp
 * \ingroup view
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::convert.
 */

#pragma once

#include <range/v3/view/transform.hpp>

#include <seqan3/core/convert.hpp>

namespace seqan3::view
{

/*!\brief Create a view that calls seqan3::convert() on each element in the input range.
 * \tparam out_t The type to convert to (must be given).
 * \param input_range The range you wish to convert, must satisfy seqan3::input_range_concept.
 * \returns A view with the value_type being `out_t`
 * \details
 * \par View properties
 * * view type: same input_range
 * * value type: out_t
 * * `const` iterable: yes
 * \par Complexity
 * Linear in the size if the input range (\f$O(n)\f$).
 * \par Exceptions
 * Strong exception guarantee (does not modify data).
 * \par Thread safety
 * Does not modify data.
 * \par Example
 * ```cpp
 * dna4_vector vec = "ACTTTGATA"_dna4;
 * auto v = vec | view::convert<char>;
 * std::cout << v << '\n'; // [A,C,T,T,T,G,A,T,A]
 * auto v2 = vec | view::convert<underlying_rank_t<dna4>> | view::convert<unsigned>;
 * std::cout << v2 << '\n'; // [0,1,3,3,3,2,0,3,0]
 *
 * std::vector<illumina18> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
 * auto v3 = qvec | view::convert<char>;
 * std::cout << v3 << '\n'; // [!,(,&,$,(,%,?,1,8]
 *
 * std::vector<dna4q> qcvec{{dna4::C, 0}, {dna4::A, 7}, {dna4::G, 5}, {dna4::T, 3}, {dna4::G, 7}, {dna4::A, 4}, {dna4::C, 30}, {dna4::T, 16}, {dna4::A, 23}};
 * auto v4 = qcvec | view::convert<char>;
 * std::cout << v4 << '\n'; // [C,A,G,T,G,A,C,T,A]
 * auto v5 = qcvec | view::convert<underlying_rank_t<dna4q>> | view::convert<unsigned>;
 * std::cout << v5 << '\n'; // [1,28,22,15,30,16,121,67,92]
 * auto v6 = qcvec | view::convert<dna4> | view::convert<char>;
 * std::cout << v6 << '\n'; // [C,A,G,T,G,A,C,T,A]
 * auto v7 = qcvec | view::convert<dna4> | view::convert<underlying_rank_t<dna4>> | view::convert<unsigned>;
 * std::cout << v7 << '\n'; // [1,0,2,3,2,0,1,3,0]
 * auto v8 = qcvec | view::convert<illumina18> | view::convert<char>;
 * std::cout << v8 << '\n'; // [!,(,&,$,(,%,?,1,8]
 * auto v9 = qcvec | view::convert<illumina18> | view::convert<underlying_rank_t<illumina18>> | view::convert<unsigned>;
 * std::cout << v9 << '\n'; // [0,7,5,3,7,4,30,16,23]
 * ```
 */
template <typename out_t>
auto const convert = ranges::view::transform([] (auto const & in) -> out_t { return seqan3::convert<out_t>(in); });

} // namespace seqan3::view

#ifndef NDEBUG
#include #include <seqan3/range/view/concept.hpp>
static_assert(seqan3::view_concept<seqan3::view::convert<char>);
#endif
