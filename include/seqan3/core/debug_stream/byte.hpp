// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
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
/*!\brief A std::byte can be printed by printing its value as integer.
 * \ingroup core_debug_stream
 */
template <>
struct std_byte_printer<std::byte>
{
    /*!\brief Prints the byte as uint8_t value.
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The output stream.
     * \param[in] arg The byte argument to print.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, std::byte const arg) const
    {
        stream << std::to_integer<uint8_t>(arg);
    }
};

} // namespace seqan3
