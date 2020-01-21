// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
 * \brief Provides the various quality score types.
 * \ingroup alphabet
 *
 * \details
 *
 * ### Introduction
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
 * ###Encoding Schemes
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
 * ###Concept
 *
 * The quality submodule defines the seqan3::writable_quality_alphabet which encompasses
 * all the alphabets, defined in the submodule, and refines the
 * seqan3::writable_alphabet by providing Phred score assignment and conversion
 * operations.
 * Additionally, this submodule defines the seqan3::quality_alphabet, which only requires
 * readablity and not assignability.
 *
 * ###Assignment and Conversion
 *
 * Quality alphabets can be converted to their char and rank representation via
 * `seqan3::to_char` and `seqan3::to_rank` respectively (like all other alphabets). Additionally
 * they can be converted to their Phred representation via `seqan3::to_phred`.
 *
 * Likewise, assignment happens via `seqan3::assign_char_to`, `seqan3::assign_rank_to` and
 * `seqan3::assign_phred_to`. Phred values outside the representable range, but inside the
 * legal range, are converted to the closest Phred score, e.g. assigning 60 to a
 * `seqan3::phred42` will result in a Phred score of 41. Assigning Phred values
 * outside the legal range results in undefined behaviour.
 *
 * All quality alphabets are explicitly convertible to each other via their
 * Phred representation. Values not present in one alphabet are mapped to the
 * closest value in the target alphabet (e.g. a `seqan3::phred63` letter with
 * value 60 will convert to a `seqan3::phred42` letter of score 41).
 */
