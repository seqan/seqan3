// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Meta-header for the \link alphabet_nucleotide Alphabet / Nucleotide submodule \endlink.
 */

/*!\defgroup alphabet_nucleotide Nucleotide
 * \brief Provides the different DNA and RNA alphabet types.
 * \ingroup alphabet
 * \see alphabet
 *
 * \details
 *
 * ### Introduction
 *
 * Nucleotide sequences are at the core of most bioinformatic data processing and while it is possible
 * to represent them in a regular std::string, it makes sense to have specialised data structures in most cases.
 * This sub-module offers multiple nucleotide alphabets that can be used with regular containers and ranges.
 *
 * | Letter   | Description            |                   seqan3::dna15        |                   seqan3::dna5         |                  seqan3::dna4          |                  seqan3::dna3bs        |                seqan3::rna15           |                    seqan3::rna5        |                 seqan3::rna4           |
 * |:--------:|------------------------|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|:--------------------------------------:|
 * |   A      | Adenine                |                              A         |                              A         |                              A         |                              A         |                              A         |                              A         |                              A         |
 * |   C      | Cytosine               |                              C         |                              C         |                              C         | <span style="color:LightGrey">T</span> |                              C         |                              C         |                              C         |
 * |   G      | Guanine                |                              G         |                              G         |                              G         |                              G         |                              G         |                              G         |                              G         |
 * |   T      | Thymine (DNA)          |                              T         |                              T         |                              T         |                              T         | <span style="color:LightGrey">U</span> | <span style="color:LightGrey">U</span> | <span style="color:LightGrey">U</span> |
 * |   U      | Uracil (RNA)           | <span style="color:LightGrey">T</span> | <span style="color:LightGrey">T</span> | <span style="color:LightGrey">T</span> | <span style="color:LightGrey">T</span> |                              U         |                              U         |                              U         |
 * |   M      | A *or* C               |                              M         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              M         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> |
 * |   R      | A *or* G               |                              R         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              R         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> |
 * |   W      | A *or* T               |                              W         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              W         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> |
 * |   Y      | C *or* T               |                              Y         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">C</span> | <span style="color:LightGrey">T</span> |                              Y         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">C</span> |
 * |   S      | C *or* G               |                              S         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">C</span> | <span style="color:LightGrey">T</span> |                              S         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">C</span> |
 * |   K      | G *or* T               |                              K         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">G</span> | <span style="color:LightGrey">G</span> |                              K         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">G</span> |
 * |   V      | A *or* C *or* G        |                              V         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              V         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> |
 * |   H      | A *or* C *or* T        |                              H         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              H         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> |
 * |   D      | A *or* G *or* T        |                              D         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              D         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">A</span> |
 * |   B      | C *or* G *or* T        |                              B         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">C</span> | <span style="color:LightGrey">T</span> |                              B         | <span style="color:LightGrey">N</span> | <span style="color:LightGrey">C</span> |
 * |   N      | A *or* C *or* G *or* T |                              N         |                              N         | <span style="color:LightGrey">A</span> | <span style="color:LightGrey">A</span> |                              N         |                              N         | <span style="color:LightGrey">A</span> |
 * | **Size** |                        |     15                                 |      5                                 |      4                                 |      3                                 |     15                                 |      5                                 |      4                                 |
 *
 * Keep in mind, that while we think of "the nucleotide alphabet" as consisting of four bases, there are indeed
 * more characters defined with different levels of ambiguity. Depending on your application it will make sense
 * to preserve this ambiguity or to discard it to save space and/or optimise computations.
 * SeqAn offers six distinct nucleotide alphabet types to accommodate for this.
 *
 * The specialised RNA alphabets are provided for convenience, however the DNA alphabets can handle being assigned a
 * ``'U'`` character, as well. See below for the details.
 *
 * Which alphabet to chose?
 *   1. in most cases, take seqan3::dna15 (includes all IUPAC characters)
 *   2. if you are memory constrained and sequence data is actually the main memory consumer, use seqan3::dna5
 *   3. if you use specialised algorithms that profit from a 2-bit representation, use seqan3::dna4
 *   4. if you are doing only RNA input/output, use the respective seqan3::rna4, seqan3::rna5, seqan3::rna15 type
 *   5. to actually save space from using smaller alphabets, you need a compressed container (e.g.
 *      seqan3::bitpacked_sequence)
 *   6. if you are working with bisulfite data use seqan3::dna3bs
 *
 * ###Printing and conversion to char
 *
 * As with all alphabets in SeqAn, none of the nucleotide alphabets can be directly converted to char or printed.
 * You need to explicitly call seqan3::to_char to convert to char. The only exception is seqan3::debug_stream
 * which does this conversion to char automatically.
 *
 * `T` and `U` are represented by the same rank and you cannot differentiate between them. The only difference between
 * e.g. seqan3::dna4 and seqan3::rna4 is the output when calling to_char().
 *
 * ###Assignment and conversions between nucleotide types
 *
 *   * Nucleotide types defined here are **implicitly** convertible to each other if they have the same size
 *     (e.g. seqan3::dna4 ↔ seqan3::rna4).
 *   * Other nucleotide types are **explicitly** convertible to each other through their character representation.
 *   * None of the nucleotide alphabets can be directly converted or assigned from `char`. You need to explicitly call
 *     `assign_char` or use a literal (see below).
 *   * Ranges of nucleotides can be converted to each other by using `std::views::transform`. See our
 *     \ref cookbook_convert_alphabet_range "cookbook" for an example.
 *
 * When assigning from `char` or converting from a larger nucleotide alphabet to a smaller one, *loss of information*
 * can occur since obviously some bases are not available. When converting to seqan3::dna5 or seqan3::rna5,
 * non-canonical bases
 * (letters other than A, C, G, T, U) are converted to ``'N'`` to preserve ambiguity at that position, while
 * for seqan3::dna4 and seqan3::rna4 they are converted to the first of the possibilities they represent (because
 * there is no letter ``'N'`` to represent ambiguity). See the greyed out values in the table at the top for
 * an overview of which conversions take place.
 *
 * `char` values that are none of the IUPAC symbols, e.g. 'P', are always converted to the equivalent of assigning 'N',
 * i.e. they result in 'A' for seqan3::dna4 and seqan3::rna4, and in 'N' for the other alphabets.
 * If the special char conversion of IUPAC characters to seqan3::dna4 is not your desired behavior, refer to our
 * cookbook for an example of \ref cookbook_custom_dna4_alphabet to change the conversion behavior.
 *
 * ###Literals
 *
 * To avoid writing ``dna4{}.assign_char('C')`` every time, you may instead use the literal ``'C'_dna4``.
 * All nucleotide types defined here have character literals (e.g \ref seqan3_dna4_char_literal "'A'_dna4") and also
 * string literals (e.g \ref seqan3_dna4_string_literal "\"ACGT\"_dna4") which return a vector of the respective type.
 *
 * ###Concept
 *
 * The nucleotide submodule defines seqan3::nucleotide_alphabet which encompasses all the alphabets defined in the
 * submodule and refines seqan3::alphabet. The only additional requirement is that their values can be
 * complemented, see below.
 *
 * ###Complement
 *
 * | Letter   | Description            | Complement |
 * |:--------:|------------------------|:----------:|
 * |   A      | Adenine                |     T      |
 * |   C      | Cytosine               |     G      |
 * |   G      | Guanine                |     C      |
 * |   T      | Thymine (DNA)          |     A      |
 * |   U      | Uracil (RNA)           |     A      |
 * |   M      | A *or* C               |     K      |
 * |   R      | A *or* G               |     Y      |
 * |   W      | A *or* T               |     W      |
 * |   Y      | C *or* T               |     R      |
 * |   S      | C *or* G               |     S      |
 * |   K      | G *or* T               |     M      |
 * |   V      | A *or* C *or* G        |     B      |
 * |   H      | A *or* C *or* T        |     D      |
 * |   D      | A *or* G *or* T        |     H      |
 * |   B      | C *or* G *or* T        |     V      |
 * |   N      | A *or* C *or* G *or* T |     N      |
 *
 * In the typical structure of DNA molecules (or double-stranded RNA), each nucleotide has a complement that it
 * pairs with. To generate the complement value of a nucleotide letter, you can call an implementation of
 * seqan3::nucleotide_alphabet::complement() on it.
 *
 * The only exception to this table is the seqan3::dna3bs alphabet. The complement for 'G' is defined as 'T' since 'C' and 'T'
 * are treated as the same letters. However, it is not recommended to use the complement of seqan3::dna3bs but rather
 * use the complement of another dna alphabet and afterwards transform it into seqan3::dna3bs.
 *
 * For the ambiguous letters, the complement is the (possibly also ambiguous) letter representing the variant of the
 * individual complements.
 *
 */

#pragma once

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna16sam.hpp>
#include <seqan3/alphabet/nucleotide/dna3bs.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
