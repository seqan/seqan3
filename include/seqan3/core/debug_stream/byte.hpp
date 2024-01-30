// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <seqan3/core/debug_stream/debug_stream_type.hpp>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
/*!\brief A std::byte can be printed by printing its value as integer.
 * \tparam    byte_type     The type of the input; must be equal to `std::byte`.
 * \param[in] s             The seqan3::debug_stream.
 * \param[in] arg           The std::byte.
 * \relates seqan3::debug_stream_type
 */
template <typename char_t, typename byte_type>
    requires std::same_as<std::remove_cvref_t<byte_type>, std::byte>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, byte_type && arg)
{
    s << std::to_integer<uint8_t>(arg);
    return s;
}

//!\}

} // namespace seqan3
