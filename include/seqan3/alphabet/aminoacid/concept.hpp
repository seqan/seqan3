// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::aminoacid_concept.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

/*!\interface seqan3::aminoacid_concept <>
 * \extends seqan3::alphabet_concept
 * \brief A concept that indicates whether an alphabet represents amino acids.
 * \ingroup aminoacid
 *
 * In addition to the requirements for seqan3::alphabet_concept, the amino_acid_concept expects
 * conforming alphabets to provide an enum-like interface with all possible 27 amino acids as values
 * (although some may be mapped to others if the alphabet is smaller than 27).
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 *
 */
//!\cond
template <typename type>
concept aminoacid_concept = requires (type v)
{
    requires alphabet_concept<type>;
    { std::remove_reference_t<type>::A } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::B } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::C } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::D } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::E } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::F } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::G } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::H } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::I } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::J } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::K } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::L } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::M } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::N } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::O } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::P } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::Q } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::R } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::S } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::T } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::U } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::V } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::W } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::X } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::Y } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::Z } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::TERMINATOR } -> std::remove_reference_t<type>;
    { std::remove_reference_t<type>::UNKNOWN } -> std::remove_reference_t<type>;

};
//!\endcond
} // namespace seqan3
