// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \sa https://en.wikipedia.org/wiki/Bit_manipulation
 */

#pragma once

#include <meta/meta.hpp>

#include <utility>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

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
 * \returns The next power of two of *n*. If *n* is already a power of two, returns *n*.
 *
 * \sa https://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 */
constexpr size_t next_power_of_two(size_t n)
{
    if (n == 0)
        return 0;

    n--;
    for(size_t shift = 1; !is_power_of_two(n + 1); shift <<= 1)
        n |= n >> shift;

    return n + 1;
}

} // namespace seqan3::detail
