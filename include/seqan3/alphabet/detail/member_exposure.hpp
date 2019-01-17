// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Free function/metafunction wrappers for alphabets with member functions/types.
 *
 * This shall not need be included manually, just include `alphabet/concept.hpp`.
 */

#pragma once

#include <iostream>
#include <optional>

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
    requires requires { typename alphabet_type_with_members::rank_type; }
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
    requires requires { alphabet_type_with_members::value_size; }
//!\endcond
struct alphabet_size<alphabet_type_with_members> :
    std::integral_constant<decltype(alphabet_type_with_members::value_size),
                           alphabet_type_with_members::value_size>
{};

/*!\brief Implementation of seqan3::semi_alphabet_concept::to_rank() that delegates to a member function.
 * \tparam alphabet_type Must provide a `.to_rank()` member function.
 * \param alph The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { { alph.to_rank() } -> underlying_rank_t<alphabet_type>; }
//!\endcond
constexpr underlying_rank_t<alphabet_type> to_rank(alphabet_type const alph) noexcept
{
    static_assert(noexcept(alph.to_rank()),
                  "The to_rank() free function can only forward to .to_rank() member functions that "
                  "are qualified as \"noexcept\".");
    return alph.to_rank();
}

/*!\brief Implementation of seqan3::semi_alphabet_concept::assign_rank() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_rank()` member function.
 * \param alph The alphabet letter that you wish to assign to.
 * \param rank The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { { alph.assign_rank(uint8_t{0}) } -> alphabet_type &; }
//!\endcond
constexpr alphabet_type & assign_rank(alphabet_type & alph, underlying_rank_t<alphabet_type> const rank) noexcept
{
    static_assert(noexcept(alph.assign_rank(rank)),
                  "The assign_rank() free function can only forward to .assign_rank() member functions that "
                  "are qualified as \"noexcept\".");
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
 * auto l = assign_rank(dna5{}, 1);  // l is of type dna5 and == 'C'_dna5
 * ~~~
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { { alph.assign_rank(uint8_t{0}) } -> alphabet_type &; }
//!\endcond
constexpr alphabet_type assign_rank(alphabet_type && alph, underlying_rank_t<alphabet_type> const rank) noexcept
{
    static_assert(noexcept(alph.assign_rank(rank)),
                  "The assign_rank() free function can only forward to .assign_rank() member functions that "
                  "are qualified as \"noexcept\".");
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
//!\cond
    requires requires (alphabet_type alph) { { alph.to_char() } -> underlying_char_t<alphabet_type>; }
//!\endcond
constexpr underlying_char_t<alphabet_type> to_char(alphabet_type const alph) noexcept
{
    static_assert(noexcept(alph.to_char()),
                  "The to_char() free function can only forward to .to_char() member functions that "
                  "are qualified as \"noexcept\".");
    return alph.to_char();
}

/*!\brief Implementation of seqan3::alphabet_concept::assign_char() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_char()` member function.
 * \param alph The alphabet letter that you wish to assign to.
 * \param chr The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { { alph.assign_char(char{0}) } -> alphabet_type &; }
//!\endcond
constexpr alphabet_type & assign_char(alphabet_type & alph, underlying_char_t<alphabet_type> const chr) noexcept
{
    static_assert(noexcept(alph.assign_char(chr)),
                  "The assign_char() free function can only forward to .assign_char() member functions that "
                  "are qualified as \"noexcept\".");
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
//!\cond
    requires requires (alphabet_type alph) { { alph.assign_char(char{0}) } -> alphabet_type &; }
//!\endcond
constexpr alphabet_type assign_char(alphabet_type && alph, underlying_char_t<alphabet_type> const chr) noexcept
{
    static_assert(noexcept(alph.assign_char(chr)),
                  "The assign_char() free function can only forward to .assign_char() member functions that "
                  "are qualified as \"noexcept\".");
    return std::move(alph.assign_char(chr));
}

/*!\brief Implementation of seqan3::alphabet_concept::char_is_valid_for() that delegates to a static member function.
 * \tparam alphabet_type Must provide a `char_is_valid()` static member function.
 * \param chr The `char` value you wish to check.
 * \returns `true` or `false`.
 */
template <typename alphabet_type_with_members>
//!\cond
    requires requires { { std::remove_reference_t<alphabet_type_with_members>::char_is_valid(char{0}) } -> bool; }
//!\endcond
constexpr bool char_is_valid_for(underlying_char_t<alphabet_type_with_members> const chr) noexcept
{
    return std::remove_reference_t<alphabet_type_with_members>::char_is_valid(chr);
}

/*!\brief Implementation of seqan3::alphabet_concept::assign_char_strict() that delegates to a member function.
 * \tparam alphabet_type Must provide an `.assign_char_strict()` member function.
 * \param alph The alphabet letter that you wish to assign to.
 * \param chr The `char` value you wish to assign.
 * \throws seqan3::invalid_char_assignment If seqan3::char_is_valid_for returns `false` for the given character.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type_with_members>
//!\cond
    requires requires (alphabet_type_with_members alph) { { alph.assign_char_strict(char{0}) } -> alphabet_type_with_members &; }
//!\endcond
constexpr alphabet_type_with_members & assign_char_strict(alphabet_type_with_members & alph,
                                                          underlying_char_t<alphabet_type_with_members> const chr)
{
    return alph.assign_char_strict(chr);
}

//!\overload
template <typename alphabet_type_with_members>
//!\cond
    requires requires (alphabet_type_with_members alph) { { alph.assign_char_strict(char{0}) } -> alphabet_type_with_members &; }
//!\endcond
constexpr alphabet_type_with_members assign_char_strict(alphabet_type_with_members && alph,
                                                        underlying_char_t<alphabet_type_with_members> const chr)
{
    return std::move(alph.assign_char_strict(chr));
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
    requires requires (nucleotide_type alph) { { alph.complement() } -> nucleotide_type; }
constexpr nucleotide_type complement(nucleotide_type const alph)
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
 * \relates seqan3::rna_structure_concept
 * \ingroup structure
 */
//!\{

/*!\brief Implementation of seqan3::rna_structure_concept::is_pair_open() that delegates to a member function.
 * \relates seqan3::rna_structure_concept
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
 * \relates seqan3::rna_structure_concept
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
 * \relates seqan3::rna_structure_concept
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

/*!\brief Specialisation of seqan3::max_pseudoknot_depth that delegates to structure_type::max_pseudoknot_depth.
 * \relates seqan3::rna_structure_concept
 * \tparam alphabet_type_with_pseudoknot_attribute Must provide a `static uint8_t max_pseudoknot_depth` member variable.
 */
template <typename alphabet_type_with_pseudoknot_attribute>
//!\cond
    requires requires (alphabet_type_with_pseudoknot_attribute)
    { { alphabet_type_with_pseudoknot_attribute::max_pseudoknot_depth } -> uint8_t; }
//!\endcond
struct max_pseudoknot_depth<alphabet_type_with_pseudoknot_attribute>
{
    //!\brief The forwarded maximum pseudoknot depth.
    static constexpr uint8_t value = alphabet_type_with_pseudoknot_attribute::max_pseudoknot_depth;
};

/*!\brief Implementation of seqan3::rna_structure_concept::pseudoknot_id() that delegates to a member function.
 * \relates seqan3::rna_structure_concept
 * \tparam alphabet_type_with_pseudoknot_attribute If it supports pseudoknots, it must provide a `.pseudoknot_id()`
 * member function, otherwise it can be omitted.
 * \param alph The alphabet letter which is checked for the pseudoknot id.
 * \returns The pseudoknot id, if alph represents an interaction, and no value otherwise.
 * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
 */
template<typename alphabet_type_with_pseudoknot_attribute>
constexpr std::optional<uint8_t> pseudoknot_id(alphabet_type_with_pseudoknot_attribute const alph)
//!\cond
    requires requires(alphabet_type_with_pseudoknot_attribute)
    { { alphabet_type_with_pseudoknot_attribute::max_pseudoknot_depth } -> uint8_t; }
//!\endcond
{
    if constexpr (alphabet_type_with_pseudoknot_attribute::max_pseudoknot_depth > 1)
        return alph.pseudoknot_id();
    else if (is_pair_open(alph) || is_pair_close(alph))
        return 0;
    else
        return std::nullopt;
}
//!\}

} // namespace seqan3
