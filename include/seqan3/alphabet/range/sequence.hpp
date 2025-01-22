// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Additional non-standard concepts for ranges.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{
/*!\interface seqan3::sequence <>
 * \brief The generic concept for a (biological) sequence.
 * \ingroup alphabet_range
 * \extends std::ranges::input_range
 *
 * A (biological) sequence is (at least) a std::ranges::input_range whose reference type models seqan3::alphabet.
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
concept sequence = std::ranges::input_range<rng_t> && alphabet<std::ranges::range_reference_t<rng_t>>;
//!\endcond
} // namespace seqan3
