// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Provides seqan3::gapped.
 */

#pragma once

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/composite/alphabet_variant.hpp>

namespace seqan3
{


/*!\brief Extends a given alphabet with a gap character.
 * \ingroup gap
 * \tparam alphabet_t Type of the letter, e.g. dna4; must satisfy seqan3::writable_alphabet.
 *
 * The gapped alphabet represents the variant of a given alphabet and the
 * seqan3::gap alphabet (e.g. the four letter DNA alphabet + a gap character).
 *
 * The gapped alphabet may be brace initialized from the static letter members of the underlying alphabet and the
 * seqan3::gap alphabet. Note that you cannot assign the alphabet by using letters of type `char`, but you instead have
 * to use the function assign_char() of the underlying alphabet or seqan3::gap::assign_char().
 *
 * \include test/snippet/alphabet/gap/gapped.cpp
 *
 * \sa For more details see alphabet_variant, which is the base class and more general than the gapped alphabet.
 */
template <typename alphabet_t>
//!\cond
    requires writable_alphabet<alphabet_t>
//!\endcond
using gapped = alphabet_variant<alphabet_t, gap>;

} // namespace seqan3

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// is_gapped_alphabet constexpr variable
// ---------------------------------------------------------------------------------------------------------------------

//!\brief Helper variable to determine if an alphabet is gapped [default: false].
//!\ingroup gap
template <typename t>
constexpr bool is_gapped_alphabet = false;

//!\brief Helper variable to determine if an alphabet is gapped, true for specilisations of seqan3::gapped.
//!\ingroup gap
template <typename t>
constexpr bool is_gapped_alphabet<gapped<t>> = true;

} // namespace seqan3::detail
