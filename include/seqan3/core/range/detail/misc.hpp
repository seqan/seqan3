// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides various utility functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Iterate over a range.
 * \ingroup core_range
 * \tparam rng_t Type of the range; must satisfy std::ranges::input_range.
 * \param rng The range.
 * \details
 * * `std::ranges::input_range<rng_t>`: Consumes input range.
 * * `std::ranges::forward_range<rng_t>`: NO-OP.
 */
template <std::ranges::input_range rng_t>
constexpr void consume(rng_t && rng)
{
    auto it = begin(rng);
    auto it_end = end(rng);
    while (it != it_end)
        ++it;
}

template <std::ranges::forward_range rng_t>
constexpr void consume(rng_t &&)
{}

} // namespace seqan3::detail
