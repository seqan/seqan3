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

/*!\file alphabet/nucleotide.hpp
 * \ingroup nucleotide
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Meta-header for the nucleotide submodule; includes all headers from alphabet/nucleotide/.
 */

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/nucl16.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

/*!\defgroup nucleotide
 * \brief Contains the different DNA and RNA alphabet types.
 * \ingroup alphabet
 *
 * \par Concept
 *
 * The nucleotide submodule defines seqan3::nucleotide_concept which encompasses all the alphabets defined in the
 * submodule and refines seqan3::alphabet_concept.
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
 * \par Alphabets
 *
 * The nucleotide module contains alphabet data types for nucleic acids. The types differ mainly
 * in their size and how they treat the T and U bases:
 *
 * | Alphabet       | Aliases                | Size | 'T' | 'U' | `UNKNOWN` |
 * |----------------|------------------------|------|-----|-----|-----------|
 * |   seqan3::dna4 |                        |    4 |  ☑  |  ☐  |       'A' |
 * |   seqan3::rna4 |                        |    4 |  ☐  |  ☑  |       'A' |
 * |   seqan3::dna5 |                        |    5 |  ☑  |  ☐  |       'N' |
 * |   seqan3::rna5 |                        |    5 |  ☐  |  ☑  |       'N' |
 * | seqan3::nucl16 | dna16, rna16, dna, rna |   16 |  ☑  |  ☑  |       'N' |
 *
 * The default and recommended alphabet is seqan3::nucl16. It can contain all the IUPAC symbols for nucleic acid data.
 * The other alphabets are smaller representations that are especially useful inside compressed containers (TODO link).
 *
 * In the smaller alphabets `T` and `U` are represented by the same rank and you cannot differentiate between them,
 * they are, however, converted into each other and not to the `UNKNOWN` character. i.e. you can use seqan3::dna5 to
 * represent RNA data, as well, the only difference is the output when calling to_char().
 *
 * \par Conversion
 *
 *   * seqan3::dna4 ↔ seqan3::rna4, as well as seqan3::dna5 ↔ seqan3::rna5 are directly assignable to each other.
 *   * All nucleotide alphabets are convertible to each other via convert().
 *   * All ranges of nucleotide alphabets are convertible to each other via seqan3::view::convert.
 */
