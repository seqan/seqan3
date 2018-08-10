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
 * \brief Genetic codes used for translating a triplet of nucleotides into an amino acid.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\brief Genetic codes used for translation of nucleotides into amino acids.
 *
 * \details
 * The numeric values of the enums correspond to the genbank transl_table values
 * (see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
 */
enum struct genetic_code : uint8_t
{
    CANONICAL=1,
//     VERT_MITOCHONDRIAL,
//     YEAST_MITOCHONDRIAL,
//     MOLD_MITOCHONDRIAL,
//     INVERT_MITOCHONDRIAL,
//     CILIATE,
//     FLATWORM_MITOCHONDRIAL = 9,
//     EUPLOTID,
//     PROKARYOTE,
//     ALT_YEAST,
//     ASCIDIAN_MITOCHONDRIAL,
//     ALT_FLATWORM_MITOCHONDRIAL,
//     BLEPHARISMA,
//     CHLOROPHYCEAN_MITOCHONDRIAL,
//     TREMATODE_MITOCHONDRIAL = 21,
//     SCENEDESMUS_MITOCHONDRIAL,
//     THRAUSTOCHYTRIUM_MITOCHONDRIAL,
//     PTEROBRANCHIA_MITOCHONDRIAL,
//     GRACILIBACTERIA
};
}
