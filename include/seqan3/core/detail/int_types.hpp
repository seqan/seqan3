// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

/*!\file
 * \brief Contains metaprogramming utilities for integer types.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \ingroup core
 */

namespace seqan3::detail
{

//!\cond DEV
//!\brief Given a value, return the smallest unsigned integer that can hold it.
template <uint64_t value>
using min_viable_uint_t = std::conditional_t<
  value <= 1ull,
  bool,
  std::conditional_t<
    value <= 255ull,
    uint8_t,
    std::conditional_t<value <= 65535ull, uint16_t, std::conditional_t<value <= 4294967295ull, uint32_t, uint64_t>>>>;

//!\brief Given a value, cast the value as the smallest unsigned integer that can hold it.
//!\sa seqan3::min_viable_uint_t
template <uint64_t value>
constexpr auto min_viable_uint_v = static_cast<min_viable_uint_t<value>>(value);
//!\endcond

} // namespace seqan3::detail
