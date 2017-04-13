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

/*!\file alphabet/view_convert.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::view::convert
 *
 * \details
 * The convert view provides easy transformation on ranges/containers/views of
 * our alphabets. This file contains the generic versions that enable
 * conversion to character and integral representations of the alphabets.
 * It also provides a fallback that tries to cast (this is e.g. helpful for
 * seqan3::alphabet_composition).
 *
 * Alphabets that wish to offer conversion to other alphabets need to explicitly
 * provide specializations.
 *
 * ~~~~~~~~~~~~~~~{.cpp}
 *  dna4_vector vec = "ACTTTGATA"_dna4;
 *  auto v = vec | view::convert<char>;
 *  std::cout << v << '\n';  // [A,C,T,T,T,G,A,T,A]
 *  auto v2 = vec | view::convert<unsigned>;
 *  std::cout << v2 << '\n'; // [0,1,3,3,3,2,0,3,0]
 *
 *  std::vector<illumina18> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
 *  auto v3 = qvec | view::convert<char>;
 *  std::cout << v3 << '\n'; // [!,(,&,$,(,%,?,1,8]
 *
 *  std::vector<dna4q> qcvec{{dna4::C, 0}, {dna4::A, 7}, {dna4::G, 5}, {dna4::T, 3}, {dna4::G, 7}, {dna4::A, 4}, {dna4::C, 30}, {dna4::T, 16}, {dna4::A, 23}};
 *  auto v4 = qcvec | view::convert<char>;
 *  std::cout << v4 << '\n'; // [C,A,G,T,G,A,C,T,A]
 *  auto v5 = qcvec | view::convert<unsigned>;
 *  std::cout << v5 << '\n'; // [1,28,22,15,30,16,121,67,92]
 *  auto v6 = qcvec | view::convert<dna4> | view::convert<char>;
 *  std::cout << v6 << '\n'; // [C,A,G,T,G,A,C,T,A]
 *  auto v7 = qcvec | view::convert<dna4> | view::convert<unsigned>;
 *  std::cout << v7 << '\n'; // [1,0,2,3,2,0,1,3,0]
 *  auto v8 = qcvec | view::convert<illumina18> | view::convert<char>;
 *  std::cout << v8 << '\n'; // [!,(,&,$,(,%,?,1,8]
 *  auto v9 = qcvec | view::convert<illumina18> | view::convert<unsigned>;
 *  std::cout << v9 << '\n'; // [0,7,5,3,7,4,30,16,23]
 * ~~~~~~~~~~~~~~~
 */

#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3::view
{

/*!\brief A view for converting ranges of seqan3::alphabet_concept. Default implementation.
 * \tparam target_type The type you wish to convert to.
 *
 * The default implementation is a last resort, it attempts a static_cast from the source
 * alphabet to the target_typ.
 */
template <typename target_type>
auto const convert = ranges::view::transform([] (alphabet_concept elem) -> target_type
{
    return static_cast<target_type>(elem);
});

/*!\brief A view for converting ranges of seqan3::alphabet_concept. Integral specialization.
 * \tparam integral_type The type you wish to convert to; requires std::is_integral_v<target_type>.
 *
 * Calls to_integral() on each letter.
 */
template <typename integral_type>
    requires std::is_integral_v<integral_type>
auto const convert<integral_type> = ranges::view::transform([] (alphabet_concept elem) -> integral_type
{
    return to_integral(elem);
});

/*!\brief A view for converting ranges of seqan3::alphabet_concept. `char` specialization.
 *
 * Calls to_char() on each letter.
 */
template <>
auto const convert<char> = ranges::view::transform([] (alphabet_concept elem) -> char { return to_char(elem); });

} // namespace seqan::view
