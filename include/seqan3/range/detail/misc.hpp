// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various utility functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Iterate over a range (consumes single-pass input ranges).
 * \ingroup range
 * \tparam rng_t Type of the range; must satisfy std::ranges::input_range.
 * \param rng The range.
 */
template <std::ranges::input_range rng_t>
constexpr void consume(rng_t && rng)
{
    auto it = begin(rng);
    auto it_end = end(rng);
    while (it != it_end)
        ++it;
}

/*!\brief Iterate over a range (NO-OP for forward ranges).
 * \ingroup range
 * \tparam rng_t Type of the range; must satisfy std::ranges::forward_range.
 */
template <std::ranges::forward_range rng_t>
constexpr void consume(rng_t &&)
{}

} // namespace seqan3::detail
