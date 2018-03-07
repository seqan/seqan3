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
 * \ingroup nucleotide
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::nucleotide_concept.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

/*!\interface seqan3::nucleotide_concept <>
 * \extends seqan3::alphabet_concept
 * \brief A concept that indicates whether an alphabet represents nucleotides.
 * \ingroup nucleotide
 *
 * In addition to the requirements for seqan3::alphabet_concept, the nucleotide_concept introduces
 * a requirement for a complement function: seqan3::nucleotide_concept::complement.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
concept bool nucleotide_concept = requires(type v)
{
    requires alphabet_concept<type>;

    {
        complement(v)
    }
    ->type;
};
//!\endcond

/*!\name Requirements for seqan3::nucleotide_concept
 * \brief You can expect these functions on all types that implement seqan3::nucleotide_concept.
 * \{
 */
/*!\fn nucleotide_type seqan3::complement(nucleotide_type const alph)
 * \brief Returns the alphabet letter's complement value.
 * \relates seqan3::nucleotide_concept
 * \param alph The alphabet letter for whom you wish to receive the complement.
 * \returns The letter's complement, e.g. 'T' for 'A'.
 * \details
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
} // namespace seqan3
