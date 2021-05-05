// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Additional non-standard concepts for ranges.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{
/*!\interface seqan3::sequence <>
 * \brief The generic concept for a (biological) sequence.
 * \ingroup range
 * \extends std::ranges::input_range
 *
 * A (biological) sequence is (at least) an std::ranges::input_range whose reference type models seqan3::alphabet.
 *
 * For example std::vector<seqan3::dna4> is a sequence of seqan3::dna4 characters.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that model this concept are shown as "implementing this interface".
 *
 * \stableapi{Since version 3.1.}
 */
//!\cond
template <typename rng_t>
SEQAN3_CONCEPT sequence = std::ranges::input_range<rng_t> && alphabet<std::ranges::range_reference_t<rng_t>>;
//!\endcond
} // namespace seqan3
