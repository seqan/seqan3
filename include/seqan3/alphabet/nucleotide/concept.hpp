// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::NucleotideAlphabet.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/std/concepts>

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

/*!\interface seqan3::NucleotideAlphabet <>
 * \extends seqan3::Alphabet
 * \brief A concept that indicates whether an alphabet represents nucleotides.
 * \ingroup nucleotide
 *
 * In addition to the requirements for seqan3::Alphabet, the NucleotideAlphabet introduces
 * a requirement for a complement function: seqan3::NucleotideAlphabet::complement.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT NucleotideAlphabet = Alphabet<type> && requires (type v, remove_cvref_t<type> c)
{
    requires std::Same<decltype(complement(v)), decltype(c)>;
};
//!\endcond

/*!\name Requirements for seqan3::NucleotideAlphabet
 * \brief You can expect these functions on all types that implement seqan3::NucleotideAlphabet.
 * \{
 */
/*!\fn nucleotide_type seqan3::complement(nucleotide_type const alph)
 * \brief Returns the alphabet letter's complement value.
 * \relates seqan3::NucleotideAlphabet
 * \param alph The alphabet letter for whom you wish to receive the complement.
 * \returns The letter's complement, e.g. 'T' for 'A'.
 * \details
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}
} // namespace seqan3
