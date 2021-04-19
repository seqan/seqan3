// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/utility/detail/bit_manipulation.hpp>
 *             instead.
 */

#pragma once

#include <cassert>

#include <seqan3/utility/detail/bits_of.hpp>
#include <seqan3/utility/detail/to_little_endian.hpp>

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0; Please #include <seqan3/std/bit> XOR <seqan3/utility/detail/bits_of.hpp> XOR <seqan3/utility/detail/to_little_endian.hpp> instead.")

namespace seqan3::detail
{
/*!\brief How many bits has a type?
 * \ingroup core
 *
 * \tparam type_t The type to determine the number of bits.
 * \deprecated This is deprecated use seqan3::detail::bits_of.
 */
template <typename type_t>
SEQAN3_DEPRECATED_310 constexpr auto sizeof_bits = seqan3::detail::bits_of<type_t>;

/*!\brief Is this number a power of two.
 * \ingroup core
 *
 * \param[in] n The number to check.
 *
 * \returns True if *n* is a power of two. False otherwise.
 *
 * \sa https://graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
 * \deprecated Use std::has_single_bit instead.
 */
SEQAN3_DEPRECATED_310 constexpr bool is_power_of_two(size_t const n)
{
    return std::has_single_bit(n);
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
 * \deprecated Use std::bit_ceil instead.
 */
SEQAN3_DEPRECATED_310 constexpr size_t next_power_of_two(size_t n)
{
    return std::bit_ceil(n);
}

/*!\brief Returns the number of 1-bits.
 * \ingroup core
 *
 * \param[in] n An unsigned integer.
 *
 * \returns The number of 1-bits.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/popcount.cpp
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Thread-safety
 *
 * Thread safe.
 *
 * ### Complexity
 *
 * Constant.
 *
 * \deprecated Use std::popcount instead.
 */
template <std::unsigned_integral unsigned_t>
SEQAN3_DEPRECATED_310 constexpr uint8_t popcount(unsigned_t const n) noexcept
{
    return std::popcount(n);
}

/*!\brief Returns the number of leading 0-bits, starting at the most significant bit position.
 * \ingroup core
 *
 * \param[in] n An unsigned integer.
 *
 * \attention *n = 0* is a special case and is undefined behaviour.
 *
 * \returns The number of leading 0-bits.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/count_leading_zeros.cpp
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Thread-safety
 *
 * Thread safe.
 *
 * ### Complexity
 *
 * Constant.
 *
 * \deprecated Use std::countl_zero instead.
 */
template <std::unsigned_integral unsigned_t>
SEQAN3_DEPRECATED_310 constexpr uint8_t count_leading_zeros(unsigned_t const n) noexcept
{
    assert(n > 0); // n == 0 might have undefined behaviour
    return std::countl_zero(n);
}

/*!\brief Returns the number of trailing 0-bits, starting at the least significant bit position.
 * \ingroup core
 *
 * \param[in] n An unsigned integer.
 *
 * \attention *n = 0* is a special case and is undefined behaviour.
 *
 * \returns The number of trailing 0-bits.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/count_trailing_zeros.cpp
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Thread-safety
 *
 * Thread safe.
 *
 * ### Complexity
 *
 * Constant.
 *
 * \deprecated Use std::countr_zero instead.
 */
template <std::unsigned_integral unsigned_t>
SEQAN3_DEPRECATED_310 constexpr uint8_t count_trailing_zeros(unsigned_t const n) noexcept
{
    assert(n > 0); // n == 0 might have undefined behaviour
    return std::countr_zero(n);
}

/*!\brief Returns the position (0-based) of the most significant bit.
 * \ingroup core
 *
 * \param[in] n An unsigned integer.
 *
 * \attention *n = 0* is a special case and is undefined behaviour.
 *
 * \returns The position of the most significant bit.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/most_significant_bit_set.cpp
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Thread-safety
 *
 * Thread safe.
 *
 * ### Complexity
 *
 * Constant.
 *
 * \deprecated Use std::bit_width(n) - 1 instead.
 */
template <std::unsigned_integral unsigned_t>
SEQAN3_DEPRECATED_310 constexpr uint8_t most_significant_bit_set(unsigned_t const n) noexcept
{
    assert(n > 0); // n == 0 might have undefined behaviour
    return std::bit_width(n) - 1;
}

} // namespace seqan3::detail
