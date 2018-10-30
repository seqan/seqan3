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
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains functions for translating a triplet of nucleotides into an amino acid.
 */

#pragma once

#include <tuple>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>
#include <seqan3/alphabet/aminoacid/translation_details.hpp>
#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/ranges>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace seqan3
{
/*!\brief Translate one nucleotide triplet into single amino acid (single nucleotide interface).
 * \tparam nucl_type The type of input nucleotides.
 * \param[in] n1 First nucleotide in triplet.
 * \param[in] n2 Second nucleotide in triplet.
 * \param[in] n3 Third nucleotide in triplet.
 *
 * \details
 *
 * Translates single nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
*/
template <genetic_code gc = genetic_code::CANONICAL, nucleotide_concept nucl_type>
constexpr aa27 translate_triplet(nucl_type const & n1, nucl_type const & n2, nucl_type const & n3) noexcept
{
    return seqan3::detail::translation_table<nucl_type, gc>::VALUE[to_rank(n1)][to_rank(n2)][to_rank(n3)];
}

/*!\brief Translate one nucleotide triplet into single amino acid (tuple interface).
 * \tparam tuple_type Type of `input_tuple`. Usually std::tuple, but similar types like std::array
 * and seqan3::pod_tuple are also supported.
 * \param[in] input_tuple Triplet of nucleotides that should be converted to amino acid.
 * \details
 *
 * Translates std::tuple or std::array with 3 nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
*/
template <genetic_code gc = genetic_code::CANONICAL, typename tuple_type>
//!\cond
    requires std::tuple_size<tuple_type>::value == 3 &&
             nucleotide_concept<std::tuple_element_t<0, tuple_type>> &&
             nucleotide_concept<std::tuple_element_t<1, tuple_type>> &&
             nucleotide_concept<std::tuple_element_t<2, tuple_type>>
//!\endcond
constexpr aa27 translate_triplet(tuple_type const & input_tuple) noexcept
{
    return translate_triplet(std::get<0>(input_tuple), std::get<1>(input_tuple), std::get<2>(input_tuple));
}

/*!\brief Translate one nucleotide triplet into single amino acid (range interface).
 * \tparam range_type Type of input_range; must satisfy std::ranges::InputRange.
 * \param[in] input_range Range of three nucleotides that should be converted to amino acid.
 *
 * \details
 *
 * Translates range with 3 nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * Strong exception guarantee (never modifies data).
*/
template <genetic_code gc = genetic_code::CANONICAL, std::ranges::InputRange range_type>
    //!\cond
    requires nucleotide_concept<std::decay_t<reference_t<std::decay_t<range_type>>>>
    //!\endcond
constexpr aa27 translate_triplet(range_type && input_range)
{
    auto n1 = begin(input_range);
    auto n2 = ++n1;
    auto n3 = ++n2;

    assert(n1 != end(input_range));
    assert(n2 != end(input_range));
    assert(n3 != end(input_range));

    return translate_triplet(*n1, *n2, *n3);
}


/*!\brief Translate one nucleotide triplet into single amino acid (range interface, input range allows random access).
 * \tparam range_type Type of input_range; must satisfy std::ranges::RandomAccessRange.
 * \param[in] input_range Range of three nucleotides that should be converted to amino acid.
 *
 * \details
 *
 * Translates range with 3 nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * Strong exception guarantee (never modifies data).
*/
template <genetic_code gc = genetic_code::CANONICAL, std::ranges::RandomAccessRange range_type>
//!\cond
    requires nucleotide_concept<std::decay_t<reference_t<std::decay_t<range_type>>>>
//!\endcond
constexpr aa27 translate_triplet(range_type && input_range)
{
    assert(input_range.begin() != end(input_range));
    assert(input_range.begin() + 1 != end(input_range));
    assert(input_range.begin() + 2 != end(input_range));

    return translate_triplet(input_range[0], input_range[1], input_range[2]);
}

} // namespace seqan3
