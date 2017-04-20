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

/*!\file alphabet/quality/aliases.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains aliases for quality_composition.
 */

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/quality/composition.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/illumina18.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/nucleotide/nucl16.hpp>

namespace seqan3
{

//!\brief An alphabet that stores a seqan3::dna4 letter and an seqan3::illumina18 letter at each position.
using dna4q = quality_composition<dna4, illumina18>;

//!\brief An alphabet that stores a seqan3::dna5 letter and an seqan3::illumina18 letter at each position.
using dna5q = quality_composition<dna5, illumina18>;

//!\brief An alphabet that stores a seqan3::rna4 letter and an seqan3::illumina18 letter at each position.
using rna4q = quality_composition<rna4, illumina18>;

//!\brief An alphabet that stores a seqan3::rna5 letter and an seqan3::illumina18 letter at each position.
using rna5q = quality_composition<rna5, illumina18>;

//!\brief An alphabet that stores a seqan3::nucl16 letter and an seqan3::illumina18 letter at each position.
using nucl16q = quality_composition<nucl16, illumina18>;

} // namespace seqan3

#ifndef NDEBUG
static_assert(seqan3::nucleotide_concept<seqan3::dna4q>);
static_assert(sizeof(seqan3::dna4q) == sizeof(seqan3::dna4) + sizeof(seqan3::illumina18));
#endif
