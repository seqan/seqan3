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
 * \brief Provides seqan3::aminoacid_concept.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

/*!\interface seqan3::aminoacid_concept <>
 * \extends seqan3::alphabet_concept
 * \brief A concept that indicates whether an alphabet represents amino acids.
 * \ingroup aminoacid
 *
 * In addition to the requirements for seqan3::alphabet_concept, the amino_acid_concept expects
 * conforming alphabets to provide an enum-like interface with all possible 27 amino acids as values
 * (although some may be mapped to others if the alphabet is smaller than 27).
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 *
 */
//!\cond
template <typename type>
concept aminoacid_concept = requires (type v)
{
    requires alphabet_concept<type>;
    { std::remove_reference_t<type>::A } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::B } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::C } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::D } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::E } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::F } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::G } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::H } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::I } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::J } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::K } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::L } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::M } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::N } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::O } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::P } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::Q } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::R } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::S } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::T } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::U } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::V } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::W } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::X } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::Y } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::Z } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::TERMINATOR } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::UNKNOWN } -> std::remove_reference_t<type>;

};
//!\endcond
} // namespace seqan3
