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
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Meta-header that includes all headers from alphabet/quality/
 */

#pragma once

#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/phred68legacy.hpp>

/*!\defgroup quality Quality
 * \brief Contains the various quality score types.
 * \ingroup alphabet
 *
 * \par Introduction
 *
 * Quality score sequences are usually output together with the DNA (or RNA)
 * sequence by sequencing machines like the Illumina Genome Analyzer. The
 * quality score of a nucleotide is also known as Phred score and is an integer
 * score being inversely proportional to the propability \f$p\f$ that a base
 * call is **incorrect**. Which roughly means that the higher a Phred score
 * is, the higher is the probabality that the corresponding nucleotide is
 * correct for that position.
 * There exists two common variants of its computation:
 *   - Salinger format with \f$Q = -10log_{10}(p)\f$
 *   - Solexa format with \f$Q = -10log_{10}(p/(1-p))\f$
 * Thus, despite implicit conversion between different quality types is
 * supported, for very low quality levels the scores vary significantly and need
 * to be corrected by an offset before being compared. The Phred score range
 * does not fit into one digit, and is therefore mapped to ASCII characters.
 * Depending on the format and the analyzer machine generation, the mappings can
 * differ.
 * Output files storing DNA sequences and their quality scores are usually
 * stored in the **FASTQ** format indicated by the file endings <I>fastq</I>
 * or <I>fq</I>.
 * This sub-module provides multiple quality alphabets that can be used in
 * combination with regular containers and ranges.
 *
 * \par Encoding Schemes
 *
 * | Format                      | Quality Type          | Phred Score Range  | Rank Range   | ASCII Range  | Assert                    |
 * |:---------------------------:|:----------------------|:------------------:|:------------:|:------------:|:-------------------------:|
 * | Sanger, Illumina 1.8+ short | seqan3::phred42       | [0 .. 41]          | [0 .. 41]    | ['!' .. 'J'] | Phred score in [0 .. 61]  |
 * | Sanger, Illumina 1.8+ long  | seqan3::phred63       | [0 .. 62]          | [0 .. 62]    | ['!' .. '_'] | Phred score in [0 .. 62]  |
 * | Solexa, Illumina [1.0; 1.8[ | seqan3::phred68legacy | [-5 .. 62]         | [0 .. 67]    | [';' .. '~'] | Phred score in [-5 .. 62] |
 *
 * The most distributed format is the *Sanger* or <I>Illumina 1.8+</I> format.
 * Despite typical Phred scores for Illumina machines range from 0 to maximal
 * 41, it is possible that processed reads reach higher scores. If you don't
 * intend to handle Phred scores larger than 41, we recommend to use
 * seqan3::phred42 due to its more space efficient implementation.
 * For other formats, like Solexa and Illumina 1.0 to 1.7 the type
 * seqan3::phred68legacy is provided. To cover also the Solexa format, the Phred
 * score is stored as a <B>signed</B> integer starting at -5.
 * An overview of all the score formats and their encodings can be found here:
 * https://en.wikipedia.org/wiki/FASTQ_format#Encoding.
 *
 * \par Concept
 *
 * The quality submodule defines the seqan3::quality_concept which encompasses
 * all the alphabets, defined in the submodule, and refines the
 * seqan3::alphabet_concept by providing Phred score assignment and conversion
 * operations.
 *
 * \par Assignment and Conversion
 *
 * Quality alphabets can be converted to their char and rank representation via
 * `to_char` and `to_rank` respectively (like all other alphabets). Additionally
 * they can be converted to their Phred representation via `to_phred`.
 *
 * Likewise, assignment happens via `assign_char`, `assign_rank` and
 * `assign_phred`. Phred values outside the representable range, but inside the
 * legal range, are converted to the closest Phred score, e.g. assigning 60 to a
 * `seqan3::phred42` will result in a Phred score of 41. Assigning Phred values
 * outside the legal range results in undefined behaviour.
 *
 * All quality alphabets are explicitly convertible to each other via their
 * Phred representation. Values not present in one alphabet are mapped to the
 * closest value in the target alphabet (e.g. a `seqan3::phred63` letter with
 * value 60 will convert to a `seqan3::phred42` letter of score 41).
 */
