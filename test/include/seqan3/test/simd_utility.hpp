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
 * \brief Contains test utilities for seqan3::simd types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/simd/all.hpp>

//!\cond DEV
/*!\brief #SIMD_EQ checks if the sizes and the content of two given
 * seqan3::simd variables matches. It is like  #EXPECT_EQ, but for seqan3::simd
 * types.
 * \ingroup simd
 * \param  left  of type seqan3::simd
 * \param  right of type seqan3::simd
 *
 * \attention
 * This macro can handle multiple "," which is normally a limitation of macros.
 *
 * \par Example
 *
 * \include test/snippet/simd/simd_test_utility.cpp
 */
#define SIMD_EQ(...) do { \
    auto [left, right] = std::make_tuple(__VA_ARGS__); \
    static_assert(seqan3::simd_concept<decltype(left)>, "The left argument of SIMD_EQ is not a simd_type"); \
    static_assert(seqan3::simd_concept<decltype(right)>, "The right argument of SIMD_EQ is not a simd_type"); \
    static_assert(std::is_same_v<decltype(left), decltype(right)>, "The left and right argument of SIMD_EQ don't have the same type."); \
    using _simd_traits_t = seqan3::simd_traits<decltype(left)>; \
    std::vector<typename _simd_traits_t::scalar_type> left_simd(_simd_traits_t::length), right_simd(_simd_traits_t::length); \
    for (size_t i = 0; i < _simd_traits_t::length; ++i) \
    std::tie(left_simd[i], right_simd[i]) = {left[i], right[i]}; \
    EXPECT_EQ(left_simd, right_simd); \
} while (false)
//!\endcond
