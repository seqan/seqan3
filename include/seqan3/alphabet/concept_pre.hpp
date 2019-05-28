// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief seqan3::Alphabet metafunction base classes.
 *
 * Include this file, if you implement an alphabet type with free/global function
 * and metafunction interfaces.
 *
 * \attention
 *
 * Note that you need to strictly follow this include order:
 * \snippet test/snippet/alphabet/concept_pre.cpp include order
 *
 * If you include `concept.hpp` before your definitions, than your type will
 * not be resolved as satisfying seqan3::Alphabet.
 *
 * This is not `true` for custom alphabets implementing the interfaces as
 * member functions/variables/types.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::RnaStructureAlphabet
// ------------------------------------------------------------------

/*!\name Requirements for seqan3::RnaStructureAlphabet
 * \brief You can expect these functions on all types that implement seqan3::RnaStructureAlphabet.
 * \relates seqan3::RnaStructureAlphabet
 * \{
 */
/*!\brief Metafunction that indicates to what extent an alphabet can handle pseudoknots.
 * [value metafunction base template]
 * \tparam alphabet_type The alphabet type whose pseudoknot ability is queried.
 * \ingroup structure
 * \details The value is the maximum allowed depth of pseudoknots.
 * A value of 1 denotes no pseudoknots `((....))`,
 * while higher values denote the maximum allowed complexity of
 * crossing interactions, e.g. depth 2 `(({....))}` or depth 3 `({[....)}]`.
 *
 * This is the expression to retrieve the value:
 * \snippet test/snippet/alphabet/concept_pre.cpp pseudoknot value retrieval
 *
 * ###Helper variable template
 *   seqan3::max_pseudoknot_depth_v as a shorthand for `seqan3::max_pseudoknot_depth<alphabet_type>::%value`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template<typename alphabet_type>
struct max_pseudoknot_depth{};

/*!\brief The pseudoknot ability of the alphabet. [value metafunction shortcut]
 * \ingroup structure
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::max_pseudoknot_depth.
 */
template<typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { max_pseudoknot_depth<alphabet_type>::value; }
//!\endcond
constexpr uint8_t max_pseudoknot_depth_v = max_pseudoknot_depth<alphabet_type>::value;
//!\}

} // namespace seqan3
