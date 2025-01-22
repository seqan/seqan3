// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 * \brief Meta-header for the \link alphabet_quality Alphabet / Quality submodule \endlink.
 */

/*!\defgroup alphabet_quality Quality
 * \brief Provides the various quality score types.
 * \ingroup alphabet
 * \see alphabet
 *
 * \details
 *
 * ### Introduction
 *
 * Quality score sequences are usually output together with the DNA (or RNA) sequence by sequencing machines like the
 * Illumina Genome Analyzer. The quality score of a nucleotide is also known as Phred score and is an integer score
 * being inversely proportional to the propability \f$p\f$ that a base call is **incorrect**. Which roughly means that
 * the higher a Phred score is, the higher is the probabality that the corresponding nucleotide is correct for that
 * position.
 * There exists two common variants of its computation:
 *   - Sanger format with \f$Q = -10\cdot\log_{10}(p)\f$
 *   - Solexa format with \f$Q = -10\cdot\log_{10}(\frac{p}{1-p})\f$
 *
 * Thus, despite implicit conversion between different quality types is supported, for very low quality levels the
 * scores vary significantly and need to be corrected by an offset before being compared.
 * For easy handling of the Phred score in file formats and console output, it is mapped to a single ASCII character.
 * The sequencing / analyser machine, e.g. HiSeq, PacBio, will dictate which Phred format is used.
 * Output files storing DNA sequences and their quality scores are usually stored in the **FASTQ** format indicated by
 * the file extensions *fastq* or *fq*.
 * This sub-module provides multiple quality alphabets that can be used in combination with regular containers and
 * ranges.
 *
 * ### Encoding Schemes
 *
 * | Standard Use Case   | Format                      | Encoding | Alphabet Type         | Phred Score Range | Rank Range |  ASCII Range                  |
 * |:-------------------:|:---------------------------:|:--------:|:----------------------|:-----------------:|:----------:|:-----------------------------:|
 * | Sanger, Illumina    | Sanger, Illumina 1.8+       | Phred+33 | seqan3::phred42       | [0 .. 41]         | [0 .. 41]  | [33 .. 74]  <br> ['!' .. 'J'] |
 * | Sanger, Illumina    | Sanger, Illumina 1.8+       | Phred+33 | seqan3::phred63       | [0 .. 62]         | [0 .. 62]  | [33 .. 95]  <br> ['!' .. '_'] |
 * | PacBio              | Sanger, Illumina 1.8+       | Phred+33 | seqan3::phred94       | [0 .. 93]         | [0 .. 93]  | [33 .. 126] <br> ['!' .. '~'] |
 * | Solexa              | Solexa, Illumina [1.0; 1.8[ | Phred+64 | seqan3::phred68solexa | [-5 .. 62]        | [0 .. 67]  | [59 .. 126] <br> [';' .. '~'] |
 *
 * The most distributed format is the *Sanger* or *Illumina 1.8+* format.
 * Despite typical Phred scores for Illumina machines range from 0 to 41, it is possible that processed reads reach
 * higher scores. If you do not intend handling Phred scores larger than 41, we recommend using seqan3::phred42 due to
 * its more space-efficient implementation (see below). If you want to store PacBio HiFi reads, we recommend to use
 * seqan3::phred94, as these use the full range of the Phred quality scores.
 * For other formats, like Solexa and Illumina 1.0 to 1.7, the type seqan3::phred68solexa is provided. To also cover the
 * Solexa format, the Phred score is stored as a **signed** integer starting at -5.
 *
 * The following figure gives a graphical explanation of the different Alphabet Types:
 *
 * <span style="font-size:12px;font-family:Courier New;line-height:1;">
        <span style="color:DarkSlateBlue"> SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS....................................................</span>\n
        <span style="color:SeaGreen">      MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM...............................</span>\n
        <span style="color:Chocolate">     PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP</span>\n
        <span style="color:BlueViolet">    ..........................OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO</span>\n\n
        <span style="font-weight:bold"> !\"\#\$\%\&'()*+,-\./0123456789:;\<\=\>?\@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{\|}~</span>\n
                                       \|.........................\|..............\|....................\|..............................\|\n
                                           33........................59............73....................95............................126\n\n
        <span style="color:DarkSlateBlue"> 0_______________________________________40.....................................................</span>\n
        <span style="color:SeaGreen">      0_______________________________________40____________________62...............................</span>\n
        <span style="color:Chocolate">     0_______________________________________40____________________62_____________________________93</span>\n
        <span style="color:BlueViolet">    .........................-5____0________9____________________________________________________62</span>\n\n
        <span style="color:DarkSlateBlue"> S - Sanger, Illumina 1.8+ - phred42</span>\n
        <span style="color:SeaGreen">      M - Sanger, Illumina 1.8+ - phred63</span>\n
        <span style="color:Chocolate">     P - Sanger, Illumina 1.8+ - phred94 (PacBio)</span>\n
        <span style="color:BlueViolet">    O - Solexa - phred68solexa</span>\n\n
        <span style="font-style:italic">Graphic was inspired by https://en.wikipedia.org/wiki/FASTQ_format#Encoding (last access 28.01.2021).</span>
   </span>
 *
 * Quality values are usually paired together with nucleotides. Therefore, it stands to reason to combine both alphabets
 * into a new data structure. In SeqAn, this can be done with seqan3::qualified. It represents the cross product between
 * a nucleotide and a quality alphabet and is the go-to choice when compression is of interest.
 *
 * The following combinations still fit into a single byte:
 * - `seqan3::qualified<seqan3::dna4, seqan3::phred42>` (alphabet size: 4 x 42 = 168)
 * - `seqan3::qualified<seqan3::dna4, seqan3::phred63>` (alphabet size: 4 x 63 = 252)
 * - `seqan3::qualified<seqan3::dna5, seqan3::phred42>` (alphabet size: 4 x 42 = 210)
 *
 * Using `seqan3::qualified` can half the storage usage compared to storing qualities and nucleotides separately. Note
 * that any combination of `seqan3::phred94` with another alphabet will cause `seqan3::qualified` to use at least 2
 * bytes.
 * While we used DNA alphabets in this example, the same properties hold true for RNA alphabets.
 *
 * ### Concept
 *
 * The quality submodule defines the seqan3::writable_quality_alphabet which encompasses all the alphabets, defined in
 * the submodule, and refines the seqan3::writable_alphabet by providing Phred score assignment and conversion
 * operations.
 * Additionally, this submodule defines the seqan3::quality_alphabet, which only requires readablity and not
 * assignability.
 *
 * ### Assignment and Conversion
 *
 * Quality alphabets can be converted to their char and rank representation via `seqan3::to_char` and `seqan3::to_rank`
 * respectively (like all other alphabets). Additionally they can be converted to their Phred representation via
 * `seqan3::to_phred`.
 *
 * Likewise, assignment happens via `seqan3::assign_char_to`, `seqan3::assign_rank_to` and `seqan3::assign_phred_to`.
 * Phred values outside the representable range, but inside the legal range, are converted to the closest Phred score,
 * e.g. assigning 60 to a `seqan3::phred42` will result in a Phred score of 41. Assigning Phred values outside the legal
 * range results in undefined behaviour.
 *
 * All quality alphabets are explicitly convertible to each other via their Phred representation. Values not present in
 * one alphabet are mapped to the closest value in the target alphabet (e.g. a `seqan3::phred63` letter with value 60
 * will convert to a `seqan3::phred42` letter of score 41, this also applies to `seqan3::phred94`).
 */

#pragma once

#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/phred68solexa.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/alphabet/quality/phred_base.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
