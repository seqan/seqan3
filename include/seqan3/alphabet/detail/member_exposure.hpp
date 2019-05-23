// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
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
// seqan3::RnaStructureAlphabet
// ------------------------------------------------------------------

/*!\name Helpers for seqan3::RnaStructureAlphabet
 * \brief These functions and metafunctions expose member variables and types so that they satisfy
 * seqan3::RnaStructureAlphabet.
 * \relates seqan3::RnaStructureAlphabet
 * \ingroup structure
 */
//!\{

/*!\brief Implementation of seqan3::RnaStructureAlphabet::is_pair_open() that delegates to a member function.
 * \relates seqan3::RnaStructureAlphabet
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

/*!\brief Implementation of seqan3::RnaStructureAlphabet::is_pair_close() that delegates to a member function.
 * \relates seqan3::RnaStructureAlphabet
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

/*!\brief Implementation of seqan3::RnaStructureAlphabet::is_unpaired() that delegates to a member function.
 * \relates seqan3::RnaStructureAlphabet
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
 * \relates seqan3::RnaStructureAlphabet
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

/*!\brief Implementation of seqan3::RnaStructureAlphabet::pseudoknot_id() that delegates to a member function.
 * \relates seqan3::RnaStructureAlphabet
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
