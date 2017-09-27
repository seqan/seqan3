// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (chr) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (chr) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Free function/metafunction wrappers for alphabets with member functions/types.
 *
 * This shall not need be included manually, just include `alphabet/concept.hpp`.
 */

#pragma once

#include <iostream>

#include <seqan3/alphabet/concept_pre.hpp>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::semi_alphabet_concept
// ------------------------------------------------------------------

/*!\name Helpers for seqan3::semi_alphabet_concept
 * \brief These functions and metafunctions expose member variables and types so that they satisfy
 * seqan3::semi_alphabet_concept.
 * \ingroup alphabet
 * \{
 */

/*!\brief Specialisation of seqan3::underlying_rank that delegates to `typename alphabet_type::rank_type`.
 * \tparam alphabet_type Must provide a member type `rank_type`.
 */
template <typename alphabet_type_with_members>
//!\cond
    requires requires (alphabet_type_with_members alph) { typename alphabet_type_with_members::rank_type; }
//!\endcond
struct underlying_rank<alphabet_type_with_members>
{
    //!\brief The forwarded rank_type.
    using type = typename alphabet_type_with_members::rank_type;
};

/*!\brief Specialisation of seqan3::alphabet_size that delegates to `alphabet_type::value_size`.
 * \tparam alphabet_type Must provide a static member variable calles `value_size`.
 *
 * Instead of accessing this struct directly, just use seqan3::alphabet_size_v.
 */
template <typename alphabet_type_with_members>
//!\cond
    requires requires (alphabet_type_with_members alph) { alphabet_type_with_members::value_size; }
//!\endcond
struct alphabet_size<alphabet_type_with_members>
{
    //!\brief The size retrieved from the type's member.
    static constexpr underlying_rank_t<alphabet_type_with_members> value =
        alphabet_type_with_members::value_size;
};

/*!\brief Implementation of seqan3::semi_alphabet_concept::to_rank() that delegates to a member function.
 * \tparam alphabet_type Must provide a `.to_rank()` member function.
 * \param alph The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename alphabet_type>
constexpr underlying_rank_t<alphabet_type> to_rank(alphabet_type const alph)
    requires requires (alphabet_type alph) { { alph.to_rank() } -> underlying_rank_t<alphabet_type>; }
{
    return alph.to_rank();
}

/*!\brief Implementation of seqan3::semi_alphabet_concept::assign_rank() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_rank()` member function.
 * \param alph The alphabet letter that you wish to assign to.
 * \param rank The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type>
constexpr alphabet_type & assign_rank(alphabet_type & alph, underlying_rank_t<alphabet_type> const rank)
    requires requires (alphabet_type alph) { { alph.assign_rank(uint8_t{0}) } -> alphabet_type &; }
{
    return alph.assign_rank(rank);
}

/*!\brief Implementation of seqan3::semi_alphabet_concept::assign_rank() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_rank()` member function.
 * \param alph An alphabet letter temporary.
 * \param rank The `rank` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from rank:
 * ~~~{.cpp}
 * auto l = assign_rank(dna5{}, 1);  // l is of type dna5 and == dna5::C
 * ~~~
 */
template <typename alphabet_type>
constexpr alphabet_type && assign_rank(alphabet_type && alph, underlying_rank_t<alphabet_type> const rank)
    requires requires (alphabet_type alph) { { alph.assign_rank(uint8_t{0}) } -> alphabet_type &; }
{
    return std::move(alph.assign_rank(rank));
}
//!\}

// ------------------------------------------------------------------
// seqan3::alphabet_concept
// ------------------------------------------------------------------

/*!\name Helpers for seqan3::alphabet_concept
 * \brief These functions and metafunctions expose member variables and types so that they satisfy
 * seqan3::alphabet_concept.
 * \ingroup alphabet
 * \{
 */

/*!\brief Specialisation of seqan3::underlying_char that delegates to `typename alphabet_type::char_type`.
 * \tparam alphabet_type Must provide a member type `char_type`.
 */
template <typename alphabet_type_with_members>
//!\cond
    requires requires (alphabet_type_with_members alph) { typename alphabet_type_with_members::char_type; }
//!\endcond
struct underlying_char<alphabet_type_with_members>
{
    //!\brief The forwarded char_type.
    using type = typename alphabet_type_with_members::char_type;
};

/*!\brief Implementation of seqan3::alphabet_concept::to_char() that delegates to a member function.
 * \tparam alphabet_type Must provide a `.to_char()` member function.
 * \param alph The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `char`).
 */
template <typename alphabet_type>
constexpr underlying_char_t<alphabet_type> to_char(alphabet_type const alph)
    requires requires (alphabet_type alph) { { alph.to_char() } -> underlying_char_t<alphabet_type>; }
{
    return alph.to_char();
}

/*!\brief Implementation of seqan3::alphabet_concept::operator<<() that delegates to `alph.to_char()`.
 * \tparam alphabet_type Must provide a `.to_char()` member function.
 * \param os The output stream you are printing to.
 * \param alph The alphabet letter that you wish to print.
 * \returns A reference to the output stream.
 */
template <typename alphabet_type>
std::ostream & operator<<(std::ostream & os, alphabet_type const alph)
//!\cond
    requires requires (alphabet_type alph) { { alph.to_char() } -> underlying_char_t<alphabet_type>; }
//!\endcond
{
    os << alph.to_char();
    return os;
}

/*!\brief Implementation of seqan3::alphabet_concept::assign_char() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_char()` member function.
 * \param alph The alphabet letter that you wish to assign to.
 * \param chr The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type>
constexpr alphabet_type & assign_char(alphabet_type & alph, underlying_char_t<alphabet_type> const chr)
    requires requires (alphabet_type alph) { { alph.assign_char(char{0}) } -> alphabet_type &; }
{
    return alph.assign_char(chr);
}

/*!\brief Implementation of seqan3::alphabet_concept::assign_char() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_char()` member function.
 * \param alph An alphabet letter temporary.
 * \param chr The `char` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from char:
 * ~~~{.cpp}
 * auto l = assign_char(dna5{}, 'G');  // l is of type dna5
 * ~~~
 */
template <typename alphabet_type>
constexpr alphabet_type && assign_char(alphabet_type && alph, underlying_char_t<alphabet_type> const chr)
    requires requires (alphabet_type alph) { { alph.assign_char(char{0}) } -> alphabet_type &; }
{
    return std::move(alph.assign_char(chr));
}
//!\}

// ------------------------------------------------------------------
// seqan3::nucleotide_concept
// ------------------------------------------------------------------

/*!\name Helpers for seqan3::nucleotide_concept
 * \brief These functions and metafunctions expose member variables and types so that they satisfy
 * seqan3::nucleotide_concept.
 * \ingroup nucleotide
 * \{
 */

/*!\brief Implementation of seqan3::nucleotide_concept::complement() that delegates to a member function.
 * \tparam nucleotide_type Must provide a `.complement()` member function.
 * \param alph The alphabet letter for whom you wish to receive the complement.
 * \returns The letter's complement, e.g. 'T' for 'A'.
 */
template <typename nucleotide_type>
constexpr nucleotide_type complement(nucleotide_type const alph)
    requires requires (nucleotide_type alph) { { alph.complement() } -> nucleotide_type; }
{
    return alph.complement();
}
//!\}

// ------------------------------------------------------------------
// seqan3::rna_structure_concept
// ------------------------------------------------------------------

/*!\name Helpers for seqan3::rna_structure_concept
 * \brief These functions and metafunctions expose member variables and types so that they satisfy
 * seqan3::rna_structure_concept.
 * \ingroup structure
 * \{
 */

/*!\brief Implementation of seqan3::rna_structure_concept::is_pair_open() that delegates to a member function.
 * \tparam structure_type Must provide a `.is_pair_open()` member function.
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents a rightward interaction, False otherwise.
 */
template <typename structure_type>
constexpr bool is_pair_open(structure_type const alph)
requires requires (structure_type alph) { { alph.is_pair_open() } -> bool; }
{
    return alph.is_pair_open();
}

/*!\brief Implementation of seqan3::rna_structure_concept::is_pair_close() that delegates to a member function.
 * \tparam structure_type Must provide a `.is_pair_close()` member function.
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents a leftward interaction, False otherwise.
 */
template <typename structure_type>
constexpr bool is_pair_close(structure_type const alph)
requires requires (structure_type alph) { { alph.is_pair_close() } -> bool; }
{
    return alph.is_pair_close();
}

/*!\brief Implementation of seqan3::rna_structure_concept::is_unpaired() that delegates to a member function.
 * \tparam structure_type Must provide a `.is_unpaired()` member function.
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents an unpaired site, False otherwise.
 */
template <typename structure_type>
constexpr bool is_unpaired(structure_type const alph)
requires requires (structure_type alph) { { alph.is_unpaired() } -> bool; }
{
    return alph.is_unpaired();
}

/*!\brief Implementation of seqan3::rna_structure_concept::pseudoknot_support() that returns a static member variable.
 * \tparam structure_type Must provide a `static bool pseudoknot_support` member variable.
 * \param alph Any alphabet letter whose type is checked for the pseudoknot property.
 * \returns True if the type can represent pseudoknots, False otherwise.
 */
template <typename structure_type>
constexpr bool pseudoknot_support(structure_type)
requires requires (structure_type) { { structure_type::pseudoknot_support } -> bool; }
{
    return structure_type::pseudoknot_support;
}
//!\}

} // namespace seqan3
