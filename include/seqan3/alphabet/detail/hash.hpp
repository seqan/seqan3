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
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides overloads for std::hash.
 */

#pragma once

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/std/ranges>

namespace std
{
/*!\brief Struct for hashing a character.
 * \ingroup alphabet
 * \tparam alphabet_t The type of character to hash; Must model seqan3::semi_alphabet_concept.
 */
template <seqan3::semi_alphabet_concept alphabet_t>
struct hash<alphabet_t>
{
    /*!\brief Compute the hash for a character.
     * \ingroup alphabet
     * \param[in] character The character to process. Must model seqan3::semi_alphabet_concept.
     *
     * \returns size_t.
     * \sa seqan3::to_rank.
     */
    size_t operator()(alphabet_t const character) const noexcept
    {
        using seqan3::to_rank;
        return to_rank(character);
    }
};

/*!\brief Struct for hashing a range of characters.
 * \ingroup alphabet
 * \tparam urng_t The type of the range; Must model std::ranges::InputRange and the reference type of the range of the
                  range must model seqan3::semi_alphabet_concept.
 */
template <ranges::InputRange urng_t>
    //!\cond
    requires seqan3::semi_alphabet_concept<seqan3::reference_t<urng_t>>
    //!\endcond
struct hash<urng_t>
{
    /*!\brief Compute the hash for a range of characters.
     * \ingroup alphabet
     * \param[in] range The input range to process. Must model std::ranges::InputRange and the reference type of the
                        range of the range must model seqan3::semi_alphabet_concept.
     * \returns size_t.
     */
    size_t operator()(urng_t const & range) const noexcept
    {
        using alphabet_t = seqan3::value_type_t<urng_t>;
        size_t result{0};
        hash<alphabet_t> h{};
        for (auto const character : range)
        {
            result *= seqan3::alphabet_size_v<alphabet_t>;
            result += h(character);
        }
        return result;
    }
};
} // namespace std
