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
 * \brief Provides seqan3::ScoringScheme.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

/*!\interface seqan3::ScoringScheme <>
 * \brief A concept that requires that type be able to score two letters.
 * \tparam t            The type the concept check is performed on (the putative scoring scheme).
 * \tparam alphabet_t   The type of the first letter that you wish to score; must model seqan3::Alphabet.
 * \tparam alphabet2_t  The type of the second letter that you wish to score; must model seqan3::Alphabet;
 *                      defaults to `alphabet_t`.
 * \ingroup scoring
 *
 * \details
 *
 * This concept makes no assumptions about configurability or assignability of the scoring scheme, only the
 * ability to score the two letters is required.
 *
 */
/*!\name Requirements for seqan3::ScoringScheme
 * \brief You can expect these members on all types that implement seqan3::ScoringScheme.
 * \memberof seqan3::ScoringScheme
 * \{
 */
/*!\typedef     typedef IMPLEMENTATION_DEFINED score_type;
 * \brief       The type returned by seqan3::ScoringScheme::score(), usually a seqan3::Arithmetic.
 * \memberof seqan3::ScoringScheme
 *
 * \details
 * \attention This is a concept requirement, not an actual typedef (however types satisfying this concept
 * will provide an implementation).
 */
/*!\fn          score_type score(alph1_t const alph1, alph2_t const alph2);
 * \brief       Compute the score of two letters
 * \param alph1 First letter.
 * \param alph2 Second letter.
 * \memberof seqan3::ScoringScheme
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
//!\cond
template <typename t, Alphabet alphabet_t, Alphabet alphabet2_t = alphabet_t>
concept ScoringScheme = requires (t scheme,
                                                alphabet_t const alph1,
                                                alphabet2_t const alph2)
{
    { scheme.score(alph1, alph2) } -> std::CommonReference<typename std::remove_reference_t<t>::score_type>;
    { scheme.score(alphabet_t{}, alphabet2_t{}) }
        -> std::CommonReference<typename std::remove_reference_t<t>::score_type>;
};
//!\endcond

} // namespace seqan3
