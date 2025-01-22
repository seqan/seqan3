// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \sa https://en.wikipedia.org/wiki/Bit_manipulation
 */

#pragma once

#include <climits>

#include <seqan3/utility/detail/integer_traits.hpp>

namespace seqan3::detail
{

/*!\brief How many bits has a type?
 * \ingroup utility
 *
 * \tparam type_t The type to determine the number of bits.
 */
template <typename type_t>
constexpr auto bits_of = min_viable_uint_v<CHAR_BIT * sizeof(type_t)>;

} // namespace seqan3::detail
