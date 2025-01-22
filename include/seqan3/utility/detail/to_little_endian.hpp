// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <bit>
#include <concepts>

// Find correct header for byte-order conversion functions.
#if __has_include(<endian.h>) // unix GLIBC
#    include <endian.h>
#elif __has_include(<sys/endian.h>) // *BSD
#    include <sys/endian.h>
#endif // __has_include(endian.h)

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\brief Convert the byte encoding of integer values to little-endian byte order.
 * \ingroup utility
 * \tparam type The type of the value to convert; must model std::integral.
 * \param  in   The input value to convert.
 * \returns the converted value in little-endian byte-order.
 *
 * \details
 *
 * This function swaps the bytes if the host system uses big endian. In this case only 1, 2, 4, or 8 byte big
 * integral types are allowed as input. On host systems with little endian this function is a no-op and returns the
 * unchanged input value. Other systems with mixed endianness are not supported.
 */
template <std::integral type>
constexpr type to_little_endian(type const in) noexcept
{
    if constexpr (std::endian::native == std::endian::little)
    {
        return in;
    }
    else if constexpr (std::endian::native == std::endian::big)
    {
        static_assert(sizeof(type) <= 8,
                      "Can only convert the byte encoding for integral numbers with a size of up to 8 bytes.");
        static_assert(std::has_single_bit(sizeof(type)),
                      "Can only convert the byte encoding for integral numbers whose byte size is a power of two.");

        if constexpr (sizeof(type) == 2)
            return htole16(in);
        else if constexpr (sizeof(type) == 4)
            return htole32(in);
        else if constexpr (sizeof(type) == 8)
            return htole64(in);
        else
            return in; // single byte.
    }
    else
    {
        static_assert(std::endian::native == std::endian::little || std::endian::native == std::endian::big,
                      "Expected a little-endian or big-endian platform.");
    }
}

} // namespace seqan3::detail
