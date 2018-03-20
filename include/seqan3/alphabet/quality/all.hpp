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

#pragma once

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Meta-header that includes all headers from alphabet/quality/
 */

#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/phred68.hpp>

/*!\defgroup quality Quality
 * \brief Contains the various quality score types.
 * \ingroup alphabet
 *
 * \par Introduction
 *
 * Quality score sequences are usually output together with the DNA (or RNA) sequence
 * by sequencing machines like the Illumina Genome Analyzer. The quality score
 * of a nucleotide is also known as Phred score and is an integer score being
 * inversely proportional to the propability p that a base call is <B>incorrect</B>.
 * Which basically means that the higher a Phred score is, the higher the probabality that the
 * corresponding nucleotide is correct for that position. There exists mainly two
 * variants of its computation:
 *   1. Salinger format with Q = -10log_10(p)
 *   2. Solexa format with Q = -10log_10(p/(1-p))
 * Thus, for very low quality levels the scores vary significantly and need to be
 * corrected by an offset before being compared. The Phred score range is not
 * displayable with one digit, and is therefore mapped to ASCII values. Depending
 * on the format and the analyzer machine generation, the mappings may differ.
 * Output files containing DNA sequences and their quality scores are usually
 * stored in the .fastq or .fq format.
 * This sub-module offers multiple quality alphabets that can be used with regular
 * containers and ranges.
 *
 * \par Encoding Schemes
 *
 * | Format                      | Quality Type | Phred Score Range  | Rank Range   | ASCII Range  |
 * |:---------------------------:|:-------------|:------------------:|:------------:|:------------:|
 * | Sanger, Illumina 1.8+ short | phred42      | [0 .. 41]          | [0 .. 41]    | ['!' .. 'J'] |
 * | Sanger, Illumina 1.8+ long  | phred63      | [0 .. 62]          | [0 .. 62]    | ['!' .. '~'] |
 * | Solexa, Illumina [1.0; 1.8[ | phred68      | [-5 .. 62]         | [0 .. 67]    | [';' .. '~'] |
 *
 * The most distributed format is the Sanger and Illumina 1.8+ format. Despite
 * typical phred scores for Illumina machines range from 0 to maximal 41, it is
 * possible that processed reads reach higher scores. If 41 is not exceeded it is
 * recommended to use phred42 to due to its smaller space consumption compared
 * to phred63.
 * For other formats, like Solexa and Illumina 1.0 to 1.7 the type phred68 is
 * provided. To cover the Solexa format, the internal rank value is a signed integer
 * starting at -5.
 * An aligned view of typical score ranges and formats can be found here:
 * https://en.wikipedia.org/wiki/FASTQ_format#Encoding.
 *
 * \par Concept
 *
 * The quality submodule defines seqan3::quality_concept which encompasses all
 * the alphabets defined in the submodule and refines seqan3::quality_concept,
 * which requires the provision of phred score assignment and conversion.
 *
 * \par Conversions
 *
 * Assignment of phred values of type phred_type is done via the =operator.
 * Except for phred68, the expected values are unsigned integers. Internally,
 * phred scores are stored as 0-based ranks of type rank_type.
 * Value assignment is possible by giving the ASCII code (assign_char), the phred
 * score (assign_phred), or directly the rank (assign_rank). Accessing the
 * different value representation is done accordingly with to_char, to_phred, and
 * to_rank.
 *
 * \par Translations
 */
