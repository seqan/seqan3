// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Meta-header for the aminoacid submodule; includes all headers from alphabet/aminoacid/.
 */
#pragma once

#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/alphabet/aminoacid/translation.hpp>
#include <seqan3/alphabet/aminoacid/translation_details.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>

/*!\defgroup aminoacid Aminoacid
 * \brief Contains the amino acid alphabets and functionality for translation from nucleotide.
 * \ingroup alphabet
 *
 * \par Introduction
 * Amino acid sequences are an important part of bioinformatic data processing and used by many applications
 * and while it is possible to represent them in a regular std::string, it makes sense to have specialised data
 * structures in most cases. This sub-module offers the 27 letter aminoacid alphabet as well as a reduced version
 * that can be used with regular container and ranges.
 * The 27 letter amino acid alphabet contains the 20 canonical amino acids, 2 additional proteinogenic amino acids
 * (Pyrrolysine and Selenocysteine) and a termination letter (*). Additionally 4 wildcard letters are offered which
 * allow a more generic usage for example in case of ambiguous amino acids (e.g. J which means either Isoleucine or
 * Leucine). See also https://en.wikipedia.org/wiki/Amino_acid for more information about the amino acid alphabet.
 *
 * \par Conversions
 * | Amino acid name            | Three letter code | One letter code | Remapped in\n seqan3::aa20      |
 * |----------------------------|-------------------|-----------------|---------------------------------|
 * |    Alanine                 | Ala               | A               | A                               |
 * |    Arginine                | Arg               | R               | R                               |
 * |    Asparagine              | Asn               | N               | N                               |
 * |    Aspartic acid           | Asp               | D               | D                               |
 * |    Cysteine                | Cys               | C               | C                               |
 * |    Tyrosine                | Tyr               | Y               | Y                               |
 * |    Glutamic acid           | Glu               | E               | E                               |
 * |    Glutamine               | Gln               | Q               | Q                               |
 * |    Glycine                 | Gly               | G               | G                               |
 * |    Histidine               | His               | H               | H                               |
 * |    Isoleucine              | Ile               | I               | I                               |
 * |    Leucine                 | leu               | L               | L                               |
 * |    Lysine                  | Lys               | K               | K                               |
 * |    Methionine              | Met               | M               | M                               |
 * |    Phenylalanine           | Phe               | F               | F                               |
 * |    Proline                 | Pro               | P               | P                               |
 * |    Serine                  | Ser               | S               | S                               |
 * |    Threonine               | Thr               | T               | T                               |
 * |    Tryptophan              | Trp               | W               | W                               |
 * |    Valine                  | Val               | V               | V                               |
 * |    Selenocysteine          | Sec               | U               | <span style="color:red">C</span>|
 * |    Pyrrolysine             | Pyl               | O               | <span style="color:red">K</span>|
 * | Asparagine or aspartic acid| Asx               | B               | <span style="color:red">D</span>|
 * | Glutamine or glutamic acid | Glx               | Z               | <span style="color:red">E</span>|
 * |    Leucine or Isoleucine   | Xle               | J               | <span style="color:red">L</span>|
 * |    Unknown                 | Xaa               | X               | <span style="color:red">S</span>|
 * |    Stop Codon              | N/A               | *               | <span style="color:red">W</span>|
 *
 * All amino acid alphabets provide static value members (like an enum) for all amino acids in the form of the
 * one-letter representation.
 * As shown above, alphabets smaller than 27 internally represent multiple amino acids as one.\n
 * For most cases it is highly recommended to use seqan3::aa27 as seqan3::aa20 provides
 * no benefits in regard to space consumption (both need 5bits).
 * Use it only when you know you need to interface with other software of formats that only support the canonical set.
 */
