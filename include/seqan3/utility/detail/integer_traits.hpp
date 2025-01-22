// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <limits>
#include <type_traits>

#include <seqan3/core/platform.hpp>

/*!\file
 * \brief Provides metaprogramming utilities for integer types.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

namespace seqan3::detail
{

// ------------------------------------------------------------------
// min_viable_uint_t
// ------------------------------------------------------------------

//!\brief Given a value, return the smallest unsigned integer that can hold it.
//!\ingroup utility
template <uint64_t value>
using min_viable_uint_t = std::conditional_t<
    value <= 1ull,
    bool,
    std::conditional_t<value <= 255ull,
                       uint8_t,
                       std::conditional_t<value <= 65535ull,
                                          uint16_t,
                                          std::conditional_t<value <= 4'294'967'295ull, uint32_t, uint64_t>>>>;

//!\brief Given a value, cast the value as the smallest unsigned integer that can hold it.
//!\ingroup utility
//!\sa seqan3::min_viable_uint_t
template <uint64_t value>
constexpr auto min_viable_uint_v = static_cast<min_viable_uint_t<value>>(value);

// ------------------------------------------------------------------
// size_in_values_v
// ------------------------------------------------------------------

//!\brief Return the number of values an integral type can have, i.e. the difference between min and max.
//!\ingroup utility
template <typename int_t>
    requires std::is_integral_v<int_t>
constexpr size_t size_in_values_v =
    static_cast<size_t>(std::numeric_limits<int_t>::max()) - std::numeric_limits<int_t>::lowest() + 1;

} // namespace seqan3::detail
