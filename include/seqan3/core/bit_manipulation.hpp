// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \sa https://en.wikipedia.org/wiki/Bit_manipulation
 */

#pragma once

#include <meta/meta.hpp>

#include <climits>
#include <utility>

// Find correct header for byte-order conversion functions.
#if __has_include(<endian.h>) // unix GLIBC
    #include <endian.h>
#elif __has_include(<sys/endian.h>)  // *BSD
    #include <sys/endian.h>
#endif // __has_include(endian.h)

#include <seqan3/core/detail/endian.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief How many bits has a type?
 * \ingroup core
 *
 * \tparam type_t The type to determine the number of bits.
 */
template <typename type_t>
constexpr auto sizeof_bits = min_viable_uint_v<CHAR_BIT * sizeof(type_t)>;

/*!\brief Is this number a power of two.
 * \ingroup core
 *
 * \param[in] n The number to check.
 *
 * \returns True if *n* is a power of two. False otherwise.
 *
 * \sa https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
 */
constexpr bool is_power_of_two(size_t const n)
{
    return n > 0 && (n & (n-1)) == 0;
}

/*!\brief Returns \f$2^{\lceil\log_2(n)\rceil}\f$ for an \f$n\f$.
 * \ingroup core
 *
 * \param[in] n A number.
 *
 * \attention *n = 0* is a special case and returns 1.
 *
 * \returns The next power of two of *n*. If *n* is already a power of two, returns *n*.
 *
 * \sa https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 */
constexpr size_t next_power_of_two(size_t n)
{
    if (n == 0)
        return 1;

    --n;
    for (size_t shift = 1; !is_power_of_two(n + 1); shift <<= 1)
        n |= n >> shift;

    return n + 1;
}

/*!\brief Returns the position of the most significant bit (counting from right to left).
 * \ingroup core
 *
 * \param[in] n An unsigned integer.
 *
 * \attention *n = 0* is a special case and is undefined behaviour.
 *
 * \returns The position of the most significant bit.
 */
template <std::UnsignedIntegral unsigned_t>
constexpr uint8_t bit_scan_reverse(unsigned_t n)
{
    assert(n > 0); // n == 0 might have undefined behaviour
#if defined(__GNUC__)
    if constexpr (sizeof(unsigned_t) == sizeof(unsigned long long))
        return sizeof_bits<unsigned long long> - __builtin_clzll(n) - 1;
    else if constexpr (sizeof(unsigned_t) == sizeof(unsigned long))
        return sizeof_bits<unsigned long> - __builtin_clzl(n) - 1;
    else
        return sizeof_bits<unsigned> - __builtin_clz(n) - 1;
#else
    uint8_t i = 0;
    for (; n != 0; n >>= 1, ++i);
    return i - 1;
#endif
}

/*!\brief Convert the byte encoding of integer values to little-endian byte order.
 * \ingroup core
 * \tparam type The type of the value to convert; must model std::Integral.
 * \param  in   The input value to convert.
 * \returns the converted value in little-endian byte-order.
 *
 * \details
 *
 * This function swaps the bytes if the host system uses big endian. In this case only 1, 2, 4, or 8 byte big
 * integral types are allowed as input. On host systems with little endian this function is a no-op and returns the
 * unchanged input value. Other systems with mixed endianness are not supported.
 */
template <std::Integral type>
constexpr type to_little_endian(type const in) noexcept
{
    if constexpr (endian::native == endian::little)
    {
        return in;
    }
    else if constexpr (endian::native == endian::big)
    {
        static_assert(sizeof(type) <= 8,
                      "Can only convert the byte encoding for integral numbers with a size of up to 8 bytes.");
        static_assert(is_power_of_two(sizeof(type)),
                      "Can only convert the byte encoding for integral numbers whose byte size is a power of two.");

        if constexpr (sizeof(type) == 2)
            return htole16(in);
        else if constexpr (sizeof(type) == 4)
            return htole32(in);
        else if constexpr (sizeof(type) == 8)
            return htole64(in);
        else
            return in;  // single byte.
    }
    else
    {
        static_assert(endian::native == endian::little || endian::native == endian::big,
                      "Expected a little-endian or big-endian platform.");
    }
}

} // namespace seqan3::detail
