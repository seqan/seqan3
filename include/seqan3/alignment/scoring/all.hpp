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
 * \brief Meta-header for the \link scoring Scoring sub-module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

/*!\defgroup scoring Scoring
 * \brief Provides the data structures used for scoring alphabets and sequences.
 * \ingroup alignment
 *
 * \details
 *
 * ### Nomenclature
 *
 * Throughout SeqAn3 we refer to the "score" of two letters or two sequences and imply that a higher score indicates
 * higher similarity and/or a closer relatedness. A lower or even negative score implies distance.
 * If you are used to dealing with "penalties" or "distances", instead think of "negative scores" when using SeqAn3
 * interfaces.
 *
 * ### Scoring two letters
 *
 * Scoring two letters of a single alphabet (or two similar alphabets) is performed by scoring schemes. A scoring
 * scheme is any type that models seqan3::scoring_scheme_concept, i.e. it must provide a member function that
 * takes the two letters and returns the scheme-specific score. Algorithms that expect a scoring scheme should check
 * this concept with their respective alphabet(s).
 *
 * Two generic scoring schemes are provided:
 *
 *   1. seqan3::nucleotide_scoring_scheme that accepts all nucleotides (and any alphabet that is explicitly
 * convertible to seqan3::dna15)
 *   2. seqan3::aminoacid_scoring_scheme that accepts all amino acids (and any alphabet that is explicitly convertible
 * to seqan3::aa27).
 *
 * These also support scoring two nucleotides/amino acids of different types and they also support modification of
 * their scores via `set_()` functions and by returning references to their internal score matrix. You can however
 * add completely different types, as long as they model seqan3::scoring_scheme_concept.
 *
 * <small><i>The base type seqan3::scoring_scheme_base is only implementation detail and not required for most users.
 * Types that model seqan3::scoring_scheme_concept can (but don't need to!) inherit from it.</i></small>
 *
 * ### Scoring gaps
 *
 * Gaps are usually not scored on a per-character-basis, because it is widely recognised that the likelihood of
 * `n` consecutive gaps is much higher than that of `n` individual gaps. This is based on the assumption that
 * one biological event often introduces more than one gap at a time and that single character gaps are not common
 * due to other bioligical factors like frame preservation in protein-coding sequences.
 *
 * Therefore the default scoring schemes do not cover the seqan3::gap alphabet and there is an additional
 * seqan3::gap_scheme whose score() function takes the length as parameter. This forces you to handle gaps separately
 * in your algorithms.
 */

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
