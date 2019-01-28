// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::rna_structure_concept.
 */

#pragma once

#include <optional>
#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{
/*!\interface seqan3::rna_structure_concept
 * \brief A concept that indicates whether an alphabet represents RNA structure.
 * \implements alphabet_concept
 * \tparam structure_type The structure alphabet type.
 * \ingroup structure
 * \details RNA structure alphabets are required to represent interactions among RNA nucleotides.
   Therefore, each structure letter can be categorised as unpaired, opening an interaction, or closing an interaction.
   Additionally, the ability of representing pseudoknots is a property of RNA structure types.
 */
/*!\fn bool is_pair_open(structure_type const alph)
 * \brief Check whether the given character represents a rightward interaction in an RNA structure.
 * \relates seqan3::rna_structure_concept
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents a rightward interaction, False otherwise.
 */
/*!\fn bool is_pair_close(structure_type const alph)
 * \brief Check whether the given character represents a leftward interaction in an RNA structure.
 * \relates seqan3::rna_structure_concept
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents a leftward interaction, False otherwise.
 */
/*!\fn bool is_unpaired(structure_type const alph)
 * \brief Check whether the given character represents an unpaired position in an RNA structure.
 * \relates seqan3::rna_structure_concept
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents an unpaired site, False otherwise.
 */
/*!\fn std::optional<uint8_t> pseudoknot_id(structure_type const alph)
 * \brief Get an identifier for a pseudoknotted interaction.
 * \relates seqan3::rna_structure_concept
 * \param alph The alphabet letter which is checked for the pseudoknot id.
 * \returns The pseudoknot id, if alph represents an interaction, and no value otherwise.
 * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
 */
/*!\struct max_pseudoknot_depth<structure_type>
 * \brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions, up to a certain depth.
 * \relates seqan3::rna_structure_concept
 * \details It is the number of distinct pairs of interaction symbols the format supports. The value 1 denotes no
 * pseudoknot support; any higher number gives the maximum nestedness. Value 0 is not allowed.
 */
//!\cond
template <typename structure_type>
SEQAN3_CONCEPT rna_structure_concept = requires(structure_type val)
{
    // requires fulfillment of alphabet concept
    requires alphabet_concept<structure_type>;

    // these are delegated to member functions, see file ../detail/member_exposure.hpp
    { is_pair_open(val) } -> bool;
    { is_pair_close(val) } -> bool;
    { is_unpaired(val) } -> bool;
    { pseudoknot_id(val) } -> std::optional<uint8_t>;

    // this is delegated to a static class variable, which must not be 0
    requires max_pseudoknot_depth<std::remove_reference_t<structure_type>>::value > 0;
    requires max_pseudoknot_depth_v<std::remove_reference_t<structure_type>> > 0;
};
//!\endcond

} // namespace seqan3
