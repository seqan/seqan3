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
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::RnaStructure.
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
/*!\interface seqan3::RnaStructure
 * \brief A concept that indicates whether an alphabet represents RNA structure.
 * \implements Alphabet
 * \tparam structure_type The structure alphabet type.
 * \ingroup structure
 * \details RNA structure alphabets are required to represent interactions among RNA nucleotides.
   Therefore, each structure letter can be categorised as unpaired, opening an interaction, or closing an interaction.
   Additionally, the ability of representing pseudoknots is a property of RNA structure types.
 */
/*!\fn bool is_pair_open(structure_type const alph)
 * \brief Check whether the given character represents a rightward interaction in an RNA structure.
 * \relates seqan3::RnaStructure
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents a rightward interaction, False otherwise.
 */
/*!\fn bool is_pair_close(structure_type const alph)
 * \brief Check whether the given character represents a leftward interaction in an RNA structure.
 * \relates seqan3::RnaStructure
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents a leftward interaction, False otherwise.
 */
/*!\fn bool is_unpaired(structure_type const alph)
 * \brief Check whether the given character represents an unpaired position in an RNA structure.
 * \relates seqan3::RnaStructure
 * \param alph The alphabet letter which is checked for the pairing property.
 * \returns True if the letter represents an unpaired site, False otherwise.
 */
/*!\fn std::optional<uint8_t> pseudoknot_id(structure_type const alph)
 * \brief Get an identifier for a pseudoknotted interaction.
 * \relates seqan3::RnaStructure
 * \param alph The alphabet letter which is checked for the pseudoknot id.
 * \returns The pseudoknot id, if alph represents an interaction, and no value otherwise.
 * It is guaranteed to be smaller than seqan3::max_pseudoknot_depth.
 */
/*!\struct max_pseudoknot_depth<structure_type>
 * \brief The ability of this alphabet to represent pseudoknots, i.e. crossing interactions, up to a certain depth.
 * \relates seqan3::RnaStructure
 * \details It is the number of distinct pairs of interaction symbols the format supports. The value 1 denotes no
 * pseudoknot support; any higher number gives the maximum nestedness. Value 0 is not allowed.
 */
//!\cond
template <typename structure_type>
concept RnaStructure = requires(structure_type val)
{
    // requires fulfillment of alphabet concept
    requires Alphabet<structure_type>;

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
