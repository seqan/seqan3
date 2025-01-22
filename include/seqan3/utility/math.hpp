// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides math related functionality.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <bit>
#include <cassert>
#include <cmath>
#include <concepts>
#include <stdexcept>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Computes the floor of the logarithm to the base of two for unsigned integers.
 * \ingroup utility
 * \param[in] n An unsigned integer.
 * \attention *n = 0* is a special case and is undefined.
 * \returns \f$ \lfloor log_2(n) \rfloor \f$.
 *
 * \details
 *
 * The difference to `std::floor(std::log2(n))` is that everything is computed *exactly* (without precision loss due to
 * promoting to `double`)
 *
 * ### Example
 *
 * \include test/snippet/utility/detail/floor_log2.cpp
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
 */
template <std::unsigned_integral unsigned_t>
constexpr unsigned_t floor_log2(unsigned_t const n) noexcept
{
    assert(n > 0u); // n == 0 is undefined behaviour
    return std::bit_width(n) - 1;
}

/*!\brief Computes the ceil of the logarithm to the base of two for unsigned integers.
 * \ingroup utility
 * \param[in] n An unsigned integer.
 * \attention *n = 0* is a special case and is undefined.
 * \returns \f$ \lceil log_2(n) \rceil \f$.
 *
 * \details
 *
 * The difference to `std::ceil(std::log2(n))` is that everything is computed *exactly* (without precision loss due to
 * promoting to `double`)
 *
 * ### Example
 *
 * \include test/snippet/utility/detail/ceil_log2.cpp
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
 */
template <std::unsigned_integral unsigned_t>
constexpr unsigned_t ceil_log2(unsigned_t const n) noexcept
{
    assert(n > 0u); // n == 0 is undefined behaviour
    return (n == 1u) ? 0u : seqan3::detail::floor_log2(n - 1u) + 1u;
}

} // namespace seqan3::detail

namespace seqan3
{
/*!\brief Computes the value of `base` raised to the power `exp`.
 * \ingroup utility
 * \param[in] base The base to compute the power for.
 * \param[in] exp The power to raise `base` to.
 * \returns \f$ base^{exp} \f$.
 * \throws std::overflow_error if an overflow occurs <b>(Only in Debug build)</b>.
 * \throws std::underflow_error if an underflow occurs <b>(Only in Debug build)</b>.
 * \sa https://en.cppreference.com/w/cpp/numeric/math/pow
 *
 * \details
 *
 * The difference to `std::pow` is that the powers of an integer `base` are computed *exact* (without precision loss due
 * to promoting to `double`) iff `exp_t` models `std::unsigned_integral` and
 *    * `base_t` models `std::unsigned_integral` (returns `uint64_t`)
 *    * `base_t` models `std::integral`, but **not** `std::unsigned_integral` (returns `int64_t`)
 *
 * In all other cases the return value and type is equivalent to that of `std::pow`.
 *
 * ### Example
 *
 * \include test/snippet/utility/pow.cpp
 */
template <typename base_t, std::unsigned_integral exp_t>
    requires (std::same_as<base_t, uint64_t> || std::same_as<base_t, int64_t>)
base_t pow(base_t base, exp_t exp)
{
    base_t result{1};
#ifndef NDEBUG
    if (base == 0)
        return 0;

    auto check = [base](base_t result)
    {
        if (base > 0)
        {
            return result > std::numeric_limits<base_t>::max() / base;
        }
        else if (result < 0) // and base < 0
        {
            return result < std::numeric_limits<base_t>::max() / base;
        }
        else // base < 0 and result > 0
        {
            return base < std::numeric_limits<base_t>::min() / result;
        }
    };

    for (exp_t i = 0; i < exp; ++i)
    {
        if (check(result))
        {
            std::string error_message{"Calculating " + std::to_string(base) + '^' + std::to_string(exp)
                                      + " will result in an "
                                      + (std::same_as<base_t, int64_t> ? "int64_t" : "uint64_t")};
            if (base < 0)
                throw std::underflow_error{error_message + " underflow."};
            else
                throw std::overflow_error{error_message + " overflow."};
        }
        result *= base;
    }
#else
    for (; exp; exp >>= 1, base *= base)
        result *= (exp & 1) ? base : 1;
#endif
    return result;
}

//!\cond
// If base and exponent are unsigned integrals, promote the base to `uint64_t`.
template <std::integral base_t, std::unsigned_integral exp_t>
    requires (std::unsigned_integral<base_t> && !std::same_as<base_t, uint64_t>)
uint64_t pow(base_t base, exp_t exp)
{
    return pow(static_cast<uint64_t>(base), exp);
}

// If the base is a signed integral and the exponent is an unsigned integral, promote the base to `int64_t`.
template <std::integral base_t, std::unsigned_integral exp_t>
    requires (!std::unsigned_integral<base_t> && !std::same_as<base_t, int64_t>)
int64_t pow(base_t base, exp_t exp)
{
    return pow(static_cast<int64_t>(base), exp);
}

// Otherwise delegate to `std::pow`.
template <typename base_t, typename exp_t>
    requires (!(std::integral<base_t> && std::unsigned_integral<exp_t>))
auto pow(base_t base, exp_t exp)
{
    return std::pow(base, exp);
}
//!\endcond

} // namespace seqan3
