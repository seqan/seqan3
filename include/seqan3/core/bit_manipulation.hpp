// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \sa https://en.wikipedia.org/wiki/Bit_manipulation
 */

#pragma once

#include <meta/meta.hpp>

#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

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

} // namespace seqan3::detail
