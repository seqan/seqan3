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

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Meta-header for the mask submodule; includes all headers from alphabet/mask/.
 */

#include <seqan3/alphabet/mask/mask.hpp>
#include <seqan3/alphabet/mask/masked_composition.hpp>

/*!\defgroup mask Mask
 * \brief Contains the mask alphabet and functionality for creating masked compositions.
 * \ingroup alphabet
 *
 * \par Introduction
 * Masks are useful as cartesian compositions when one wants to create a masked alphabet with
 * don't care positions, but does not want to use the seqan3::dna15::N or
 * seqan3::aa27::X because of loss of information. It will instead mark the specified characters as masked,
 * and display them as lowercase representations when printed.\n
 * There are two types of masking: "hard-masking" which converts to the UNKNOWN character and
 * "soft-masking", which is visualised by using lower-case instead of upper-case.
 * However because regular nucleotide and aminoacid alphabets discard case on assignment,
 * one needs to create additional alphabets to preserve this information (if desired).\n
 * This alphabet in itself is not useful to users directly, but instead the composition seqan3::masked may be used to
 * transform another alphabet into a new alphabet that can represent the original alphabet plus masking information.
 */
