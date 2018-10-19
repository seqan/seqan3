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

#if __has_include(<umesimd/UMESimd.h>)

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/ume_simd.hpp>

#include <iostream>

using namespace seqan3;
using namespace seqan3::detail;

using int16x8_t = UME::SIMD::SIMDVec<int16_t, 8u>;
using int32x4_t = UME::SIMD::SIMDVec<int32_t, 4u>;
using int64x2_t = UME::SIMD::SIMDVec<int64_t, 2u>;

using uint16x16_t = UME::SIMD::SIMDVec<uint16_t, 16u>;
using uint32x8_t = UME::SIMD::SIMDVec<uint32_t, 8u>;
using uint64x4_t = UME::SIMD::SIMDVec<uint64_t, 4u>;

using mask2_t  = UME::SIMD::SIMDVecMask<2>;
using mask4_t  = UME::SIMD::SIMDVecMask<4>;
using mask8_t  = UME::SIMD::SIMDVecMask<8>;
using mask16_t  = UME::SIMD::SIMDVecMask<16>;

using swizzle2_t  = UME::SIMD::SIMDSwizzle<2>;
using swizzle4_t  = UME::SIMD::SIMDSwizzle<4>;
using swizzle8_t  = UME::SIMD::SIMDSwizzle<8>;
using swizzle16_t  = UME::SIMD::SIMDSwizzle<16>;

TEST(ume_simd, ume_simd)
{
    EXPECT_TRUE((std::is_same_v<ume_simd<int16_t, 8>::type, int16x8_t>));
    EXPECT_TRUE((std::is_same_v<ume_simd<int32_t, 4>::type, int32x4_t>));
    EXPECT_TRUE((std::is_same_v<ume_simd<int64_t, 2>::type, int64x2_t>));

    EXPECT_TRUE((std::is_same_v<ume_simd<uint16_t, 16>::type, uint16x16_t>));
    EXPECT_TRUE((std::is_same_v<ume_simd<uint32_t, 8>::type, uint32x8_t>));
    EXPECT_TRUE((std::is_same_v<ume_simd<uint64_t, 4>::type, uint64x4_t>));
}

TEST(ume_simd, is_ume_simd)
{
    EXPECT_FALSE(is_ume_simd<short>::value);
    EXPECT_FALSE(is_ume_simd<int>::value);
    EXPECT_FALSE(is_ume_simd<int[15]>::value);
    EXPECT_FALSE(is_ume_simd<int*>::value);

    EXPECT_TRUE(is_ume_simd<int16x8_t>::value);
    EXPECT_TRUE(is_ume_simd<int32x4_t>::value);
    EXPECT_TRUE(is_ume_simd<int64x2_t>::value);

    EXPECT_TRUE(is_ume_simd<uint16x16_t>::value);
    EXPECT_TRUE(is_ume_simd<uint32x8_t>::value);
    EXPECT_TRUE(is_ume_simd<uint64x4_t>::value);
}

TEST(ume_simd, simd_traits)
{
    EXPECT_TRUE((std::is_same_v<simd_traits<int16x8_t>::scalar_type, int16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int32x4_t>::scalar_type, int32_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int64x2_t>::scalar_type, int64_t>));

    EXPECT_EQ(simd_traits<int16x8_t>::length, 8u);
    EXPECT_EQ(simd_traits<int32x4_t>::length, 4u);
    EXPECT_EQ(simd_traits<int64x2_t>::length, 2u);

    EXPECT_EQ(simd_traits<int16x8_t>::max_length, 16u);
    EXPECT_EQ(simd_traits<int32x4_t>::max_length, 16u);
    EXPECT_EQ(simd_traits<int64x2_t>::max_length, 16u);

    EXPECT_TRUE((std::is_same_v<simd_traits<int16x8_t>::mask_type, mask8_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int32x4_t>::mask_type, mask4_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int64x2_t>::mask_type, mask2_t>));

    EXPECT_TRUE((std::is_same_v<simd_traits<int16x8_t>::swizzle_type, swizzle8_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int32x4_t>::swizzle_type, swizzle4_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int64x2_t>::swizzle_type, swizzle2_t>));

    // avx2 (256bit)

    EXPECT_TRUE((std::is_same_v<simd_traits<uint16x16_t>::scalar_type, uint16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint32x8_t>::scalar_type, uint32_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint64x4_t>::scalar_type, uint64_t>));

    EXPECT_EQ(simd_traits<uint16x16_t>::length, 16u);
    EXPECT_EQ(simd_traits<uint32x8_t>::length, 8u);
    EXPECT_EQ(simd_traits<uint64x4_t>::length, 4u);

    EXPECT_EQ(simd_traits<uint16x16_t>::max_length, 32u);
    EXPECT_EQ(simd_traits<uint32x8_t>::max_length, 32u);
    EXPECT_EQ(simd_traits<uint64x4_t>::max_length, 32u);

    EXPECT_TRUE((std::is_same_v<simd_traits<uint16x16_t>::mask_type, mask16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint32x8_t>::mask_type, mask8_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint64x4_t>::mask_type, mask4_t>));

    EXPECT_TRUE((std::is_same_v<simd_traits<uint16x16_t>::swizzle_type, swizzle16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint32x8_t>::swizzle_type, swizzle8_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint64x4_t>::swizzle_type, swizzle4_t>));
}

TEST(ume_simd, default_simd_max_length)
{
    constexpr auto default_simd_max_length_v = default_simd_max_length<ume_simd>;
#if defined(__AVX512F__)
    EXPECT_EQ(default_simd_max_length_v, 64u);
#elif defined(__AVX2__)
    EXPECT_EQ(default_simd_max_length_v, 32u);
#elif defined(__SSE4_2__)
    EXPECT_EQ(default_simd_max_length_v, 16u);
#else
    EXPECT_EQ(default_simd_max_length_v, 0u);
#endif
}

TEST(ume_simd, simd_concept)
{
    EXPECT_FALSE(simd_concept<short>);
    EXPECT_FALSE(simd_concept<int>);
    EXPECT_FALSE(simd_concept<int[15]>);
    EXPECT_FALSE(simd_concept<int*>);

    EXPECT_TRUE(simd_concept<int16x8_t>);
    EXPECT_TRUE(simd_concept<int32x4_t>);
    EXPECT_TRUE(simd_concept<int64x2_t>);
    EXPECT_TRUE(simd_concept<uint16x16_t>);
    EXPECT_TRUE(simd_concept<uint32x8_t>);
    EXPECT_TRUE(simd_concept<uint64x4_t>);
}

#else // __has_include(<umesimd/UMESimd.h>)

TEST(ume_simd, DISABLED_not_included)
{}

#endif // __has_include(<umesimd/UMESimd.h>)
