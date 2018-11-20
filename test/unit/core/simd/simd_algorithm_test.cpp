// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <gtest/gtest.h>

#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/test/simd_utility.hpp>

#include <iostream>

using namespace seqan3;

TEST(simd_algorithm, fill)
{
    using simd_type = simd_type_t<int16_t, 8>;

    constexpr simd_type result = simd::fill<simd_type>(4);
    simd_type expect = detail::simd_transform<simd_type>([](size_t)
    {
        return 4;
    });
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, iota)
{
    using simd_type = simd_type_t<int16_t, 8>;

    constexpr simd_type result = simd::iota<simd_type>(0);
    simd_type expect = detail::simd_transform<simd_type>([i = 0](size_t) mutable
    {
        return i++;
    });
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, iota_with_offset)
{
    using simd_type = simd_type_t<int16_t, 8>;

    constexpr simd_type result = simd::iota<simd_type>(5);
    simd_type expect = detail::simd_transform<simd_type>([](size_t i)
    {
        return 5 + i;
    });
    SIMD_EQ(result, expect);
}
