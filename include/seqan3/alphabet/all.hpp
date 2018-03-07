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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Meta-header for the alphabet module.
 *
 * \defgroup alphabet Alphabet
 *
 * ## Introduction
 *
 * Alphabets are a core component in SeqAn. They enable us to represent the smallest
 * unit of biological sequence data, e.g. a nucleotide or an amino acid.
 *
 * In theory, these could just be represented as a `char` and this is how many people
 * perceive them, but it makes sense to use a smaller, stricter and well-defined alphabet
 * in almost all cases, because:
 *
 *   * Most of our sequences alphabets are actually smaller and can possibly be **represented by less bits**
 * than a `char`, e.g. a `char` can have 256 values and thus must be represented by 8 bits of memory, but a DNA
 * character *could* be represented by 2 bits, because it only has four values in the smallest representation
 * ('A', 'C', 'G', 'T').
 *   * From a programmer's point of view it is very important to also access the **rank of a letter**, i.e. we need to
 * be able to convert 'A', 'C', 'G', 'T' to `0`, `1`, `2`, `3` respectively. In fact the rank representation is used a
 * lot more often than the visual representation which is only used in input/output.
 *   * You may want to prevent the user from selecting **letters that aren't part of that alphabet**, e.g. 'J' is
 * not part of any nucleotide definition so people shouldn't be able to select 'J'. With a specialised
 * alphabet you could instead convert those to a special "unknown"-character, e.g. you might want 'J' to be
 * always converted to 'N'.
 *
 * In SeqAn there are alphabet types for typical sequence alphabets like DNA and amino acid, but also
 * for qualities, RNA structures and alignment gaps. In addition there are templates for combining alphabet
 * types into new alphabets, and wrappers for existing data types like the canonical `char`.
 *
 * To be included into the alphabet module, an alphabet must satisfy the generic seqan3::alphabet_concept
 * documented below. While this only encompasses a minimum set of requirements, many of our alphabets provide
 * more features and there are more refined concepts. The inheritance diagram of seqan3::alphabet_concept gives
 * a detailed overview. A more basic overview of this module and it's submodules is available in the collaboration
 * diagram at the top of this page.
 *
 * ## The alphabet concept
 *
 * The seqan3::alphabet_concept defines the requirements a type needs to meet to be considered an alphabet
 * by SeqAn, or in other words: you can expect certain properties and functions to be defined on
 * all data types we call an alphabet.
 *
 * ### Alphabet size
 *
 * All alphabets in SeqAn have a fixed size. It
 * can be queried via the seqan3::alphabet_size metafunction and *optionally* also the `value_size` static
 * member of the alphabet (see below for "members VS free/global functions").
 *
 * In some areas we provide alphabets types with different sizes for the same purpose, e.g. seqan3::dna4
 * ('A', 'C', 'G', 'T'), seqan3::dna5 (plus 'N') and seqan3::dna15 (plus ambiguous characters defined by
 * IUPAC). By convention most of our alphabets carry their size in their name (seqan3::dna4 has size 4 a.s.o.).
 *
 * A main reason for choosing a smaller alphabet over a bigger one is the possibility of **optimising for
 * space efficiency**. Note, however, that a single letter by itself can never be smaller than a byte for
 * architectural reasons. Actual space improvements are realised via secondary structures, e.g. when
 * using a seqan3::bitcompressed_vector<seqan3::dna4> instead of std::vector<seqan3::dna4>. Also
 * the single letter quality composition seqan3::qualified<seqan3::dna4, seqan3::phred42> fits into one byte, because
 * the product of the alphabet sizes (4 * 42) is smaller than 256; whereas the same composition
 * with seqan3::dna15 requires two bytes per letter (15 * 42 > 256).
 *
 * ### Assigning and retrieving values
 *
 * As mentioned above, we typically think of alphabets in their character representation, but we also
 * require them in "rank representation" as programmers. In C and C++ it is quite difficult to
 * cleanly differentiate between these, because the `char` type is considered an integral type and can
 * be used to index an array (e.g. <tt>my_array['A']</tt> translates to `my_array[65]`). Moreover the sign of `char` is
 * implementation defined and on many platforms the smallest integer types `int8_t` and `uint8_t` are literally the
 * same types as `signed char` and `unsigned char` respectively.
 *
 * This leads to ambiguity when assigning and retrieving values:
 *
 * ```cpp
 * // WRONG EXAMPLE:
 * dna4 my_letter{0};      // we want to set the default, an A
 * dna4 my_letter{'A'};    // we also want to set an A, but we are setting value 65
 *
 * std::cout << my_letter; // you expect 'A', but how would you access the number?
 * ```
 *
 * To solve this problem, every alphabet defines two interfaces:
 *
 *   1. a **character based interface** with
 *     1. the underlying character type able to represent this alphabet visually (almost always `char`,
 * but could be `char16_t` or `char32_t`, as well)
 *     2. a seqan3::to_char function to produce the visual representation
 *     3. a seqan3::assign_char function to assign from the visual representation
 *   2. a **rank based interface** with
 *     1. the underlying rank type able to represent this alphabet numerically; this type must be able to represent
 * the numbers from `0` to `alphabet size - 1` (often `uint8_t`, but sometimes a larger unsigned integral type)
 *     2. a seqan3::to_rank function to produce the numerical representation
 *     3. a seqan3::assign_rank function to assign from the numerical representation
 *
 * To prevent the aforementioned ambiguity, you can neither assign from rank or char representation via `operator=`,
 * nor can you cast the alphabet to either of it's representation forms, **you need to explicitly use the
 * interfaces**:
 * ```cpp
 * dna4 my_letter;
 * assign_rank(my_letter, 0);       // assign an A via rank interface
 * assign_char(my_letter, 'A');     // assign an A via char interface
 * my_letter = dna4::A;             // some alphabets(BUT NOT ALL!) also provide an enum-like interface
 *
 * std::cout << to_char(my_letter);            // prints 'A'
 * std::cout << (unsigned)to_rank(my_letter);  // prints 0
 * // we have to add the cast here, because uint8_t is also treated as a char type by default :(
 * ```
 *
 * For efficiency, the representation saved internally is normally the rank representation, and the character
 * representation
 * is generated via conversion tables. This is, however, not required as long as both interfaces are provided
 * and all functions operate in constant time.
 *
 * <small>In the documentation you will also encounter seqan3::semi_alphabet_concept. It describes "one half" of an
 * alphabet and only defines the rank interface as a type requirement. It is mainly used internally and not
 * relevant to most users of SeqAn.</small>
 *
 * ### Members VS free/global functions
 *
 * The alphabet concept (as most concepts in SeqAn) looks for free/global functions, i.e. you need to be able
 * to call `seqan3::to_rank(my_letter)`, however *most* alphabets also provide a member function, i.e.
 * `my_letter.to_rank()`. The same is true for the metafunction seqan3::alphabet_size vs the static data member
 * `value_size`.
 *
 * Members are provided for convenience and if you are an application developer who works with a single concrete
 * alphabet type you are fine with using the member functions. If you, however, implement a generic function
 * that accepts different alphabet types, you need to use the free function / metafunction interface, because
 * it is the only interface guaranteed to exist (member functions are **not** required/enforced by the concept).
 *
 * ## Containers over alphabets
 *
 * In SeqAn3 it is recommended you use the STL container classes like std::vector for storing sequence data,
 * but you can use other class templates if they satisfy the respective seqan3::container_concept, e.g. `std::deque` or
 * <a href="https://github.com/facebook/folly/blob/master/folly/docs/FBVector.md" target="_blank">
 * <tt>folly::fbvector</tt></a> or even <a href="http://doc.qt.io/qt-5/qvector.html" target="_blank">
 * <tt>Qt::QVector</tt></a>.
 *
 * `std::basic_string` is also supported, however, we recommend against using it,
 * because it is not safe (and not useful) to call certain members like `.c_str()` if our alphabets are used as
 * value type.
 *
 * We provide specialised containers with certain properties in the \ref range module.
 *
 */

#pragma once

#include <seqan3/alphabet/adaptation/all.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/composition/all.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/alphabet/structure/all.hpp>
