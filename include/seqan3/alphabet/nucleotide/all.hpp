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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Meta-header for the nucleotide submodule; includes all headers from alphabet/nucleotide/.
 */

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>

/*!\defgroup nucleotide Nucleotide
 * \brief Contains the different DNA and RNA alphabet types.
 * \ingroup alphabet
 *
 * \par Introduction
 *
 * Nucleotide sequences are at the core of most bioinformatic data processing and while it is possible
 * to represent them in a regular std::string, it makes sense to have specialised data structures in most cases.
 * This sub-module offers multiple nucleotide alphabets that can be used with regular containers and ranges.
 *
 * | Letter   | Description            | seqan3::dna15  | seqan3::dna5 | seqan3::dna4 | seqan3::rna15  | seqan3::rna5 | seqan3::rna4 |
 * |:--------:|------------------------|:--------------:|:------------:|:------------:|:--------------:|:------------:|:------------:|
 * |   A      | Adenine                |      ☑         |      ☑       |      ☑       |      ☑         |      ☑       |      ☑       |
 * |   C      | Cytosine               |      ☑         |      ☑       |      ☑       |      ☑         |      ☑       |      ☑       |
 * |   G      | Guanine                |      ☑         |      ☑       |      ☑       |      ☑         |      ☑       |      ☑       |
 * |   T      | Thymine (DNA)          |      ☑         |      ☑       |      ☑       |      ☐         |      ☐       |      ☐       |
 * |   U      | Uracil (RNA)           |      ☐         |      ☐       |      ☐       |      ☑         |      ☑       |      ☑       |
 * |   M      | A *or* C               |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   R      | A *or* G               |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   W      | A *or* T               |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   Y      | C *or* T               |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   S      | C *or* G               |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   K      | G *or* T               |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   V      | A *or* C *or* G        |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   H      | A *or* C *or* T        |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   D      | A *or* G *or* T        |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   B      | C *or* G *or* T        |      ☑         |      ☐       |      ☐       |      ☑         |      ☐       |      ☐       |
 * |   N      | A *or* C *or* G *or* T |      ☑         |      ☑       |      ☐       |      ☑         |      ☑       |      ☐       |
 * | **Size** |                        |     15         |      5       |      4       |     15         |      5       |      4       |
 *
 * Keep in mind, that while we think of "the nucleotide alphabet" as consisting of four bases, there are indeed
 * more characters defined with different levels of ambiguity. Depending on your application it will make sense
 * to preserve this ambiguity or to discard it to save space and/or optimise computations.
 * SeqAn offers six distinct nucleotide alphabet types to accommodate for this.
 *
 * The specialised RNA alphabets are provided for convenience, however the DNA alphabets provide a `::%U` member
 * and can handle being assigned a `'U'` character, as well. See below for the details.
 *
 * Which alphabet to chose?
 *   1. in most cases, take seqan3::dna15
 *   2. if you are memory constrained and sequence data is actually the main memory consumer, use seqan3::dna5
 *   3. if you use specialised algorithms that profit from a 2-bit representation, use seqan3::dna4
 *   4. if you are doing only RNA input/output, use the respective seqan3::rna* type
 *   5. to actually save space from using smaller alphabets, you need a compressed container (TODO link)
 *
 * \par Concept
 *
 * The nucleotide submodule defines seqan3::NucleotideAlphabet which encompasses all the alphabets defined in the
 * submodule and refines seqan3::Alphabet.
 *
 * While it is not required by the concept, all alphabets declared in the module offer an `enum`-like interface
 * to the individual letters so you can use them as values in assignments and comparisons:
 *
 * ~~~{.cpp}
 * dna4 l{dna4::A};
 *
 * if (l == dna4::C}
 *     // ...
 * ~~~
 *
 * All alphabets define at least the letters `::%A`, `::%C`, `::%G`, `::%T`, `::%U`, `::%UNKNOWN` (although some
 * might be aliases of others).
 *
 * \par Printing and conversion to char
 *
 * | to_char()      | seqan3::dna15  | seqan3::dna5 | seqan3::dna4 | seqan3::rna15  | seqan3::rna5 | seqan3::rna4 |
 * |:--------------:|:--------------:|:------------:|:------------:|:--------------:|:------------:|:------------:|
 * | type::A        |    `'A'`       |    `'A'`     |    `'A'`     |    `'A'`       |    `'A'`     |    `'A'`     |
 * |   ...          |    ...         |    ...       |    ...       |    ...         |    ...       |    ...       |
 * | type::T        |    `'T'`       |    `'T'`     |    `'T'`     |    `'U'`       |    `'U'`     |    `'U'`     |
 * | type::U        |    `'T'`       |    `'T'`     |    `'T'`     |    `'U'`       |    `'U'`     |    `'U'`     |
 * | type::R        |    `'R'`       |     ---      |    ---       |    `'R'`       |     ---      |     ---      |
 * |   ...          |    ...         |     ---      |    ---       |    ...         |     ---      |     ---      |
 * | type::N        |    `'N'`       |    `'N'`     |    ---       |    `'N'`       |    `'N'`     |     ---      |
 * | type::UNKNOWN  |    `'N'`       |    `'N'`     |    `'A'`     |    `'N'`       |    `'N'`     |    `'A'`     |
 *
 * `T` and `U` are represented by the same rank and you cannot differentiate between them
 * (the static members `::%T` and `::%U` are synonymous). The only difference between e.g. seqan3::dna4 and
 * seqan3::rna4 is the output when calling to_char().
 *
 * \par Assignment and conversions between nucleotide types
 *
 *   * Nucleotide types are **implicitly** convertible to each other if they have the same size
 * (e.g. seqan3::dna4 ↔ seqan3::rna4).
 *   * Other nucleotide types are **explicitly** convertible to each other through their character representation.
 *   * All ranges of nucleotide alphabets are convertible to each other via seqan3::view::convert.
 *
 * | assign_char() | seqan3::dna15  | seqan3::dna5  | seqan3::dna4  | seqan3::rna15  | seqan3::rna5 | seqan3::rna4  |
 * |:-------------:|:--------------:|:-------------:|:-------------:|:--------------:|:------------:|:-------------:|
 * |   `'A'`       |   dna15::A     |  dna5::A      | dna4::A       |   rna15::A     |  rna5::A     | rna4::A       |
 * |   ...         |    ...         |    ...        | ...           |    ...         |    ...       |     ...       |
 * |   `'T'`       |   dna15::T ¹   |  dna5::T ¹    | dna4::T ¹     |   rna15::T ¹   |  rna5::T ¹   | rna4::T ¹     |
 * |   `'U'`       |   dna15::T ¹   |  dna5::T ¹    | dna4::T ¹     |   rna15::T ¹   |  rna5::T ¹   | rna4::T ¹     |
 * |   `'R'`       |   dna15::R     |  dna5::N ²    | dna4::A ³     |   rna15::R     |  rna5::N ²   | rna4::A ³     |
 * |   `'Y'`       |   dna15::Y     |  dna5::N ²    | dna4::C ³     |   rna15::Y     |  rna5::N ²   | rna4::C ³     |
 * |   ...         |    ...         |  dna5::N ²    | ... ³         |    ...         |  rna5::N ²   | ... ³         |
 * |   `'N'`       |   dna15::N     |  dna5::N ²    | dna4::A ³     |   rna15::N     |  rna5::N ²   | rna4::A ³     |
 * | e.g. `'!'`    |   dna15::N     |  dna5::N ²    | dna4::A ²     |   rna15::N     |  rna5::N ²   | rna4::A ²     |
 *
 * ¹ same as `type::U`<br>
 * ² same as `type::UNKNOWN`<br>
 * ³ the first character in the ambiguous set
 *
 * When converting from `char` or from a larger nucleotide alphabet to a smaller one, *loss of information* can occur
 * since obviously some bases are not available. When converting to seqan3::dna5 or seqan3::rna5, non-canonical bases
 * (letters other than A, C, G, T, U) are converted to `'N'` to preserve ambiguity at that position, while
 * for seqan3::dna4 and seqan3::rna4 they are converted to the first of the possibilities they represent (because
 * there is no letter `'N'` to represent ambiguity).
 *
 * \par Complement
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
 * |   Y      | C *or* T               |     S      |
 * |   S      | C *or* G               |     R      |
 * |   K      | G *or* T               |     M      |
 * |   V      | A *or* C *or* G        |     B      |
 * |   H      | A *or* C *or* T        |     D      |
 * |   D      | A *or* G *or* T        |     H      |
 * |   B      | C *or* G *or* T        |     V      |
 * |   N      | A *or* C *or* G *or* T |     N      |
 *
 * In the typical structure of DNA molecules (or double-stranded RNA), each nucleotide has a complement that it
 * pairs with. To generate the complement value of a nucleotide letter, you can call an implementation of
 * seqan3::NucleotideAlphabet::complement() on it.
 *
 * For the ambiguous letters, the complement is the (possibly also ambiguous) letter representing the union of the
 * individual complements.
 *
 */
