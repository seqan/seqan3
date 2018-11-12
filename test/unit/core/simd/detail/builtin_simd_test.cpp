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

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/builtin_simd.hpp>

#include <iostream>
#include <type_traits>

using namespace seqan3;

#ifdef __clang__
using uint8x16_t = uint8_t __attribute__((ext_vector_type(16)));
using int16x8_t = int16_t __attribute__((ext_vector_type(8)));
using int32x4_t = int32_t __attribute__((ext_vector_type(4)));
using int64x2_t = int64_t __attribute__((ext_vector_type(2)));

using int16x16_t = int16_t __attribute__((ext_vector_type(16)));
using int32x8_t = int32_t __attribute__((ext_vector_type(8)));
using int64x4_t = int64_t __attribute__((ext_vector_type(4)));

using uint8x32_t = uint8_t __attribute__((ext_vector_type(32)));
using uint16x16_t = uint16_t __attribute__((ext_vector_type(16)));
using uint32x8_t = uint32_t __attribute__((ext_vector_type(8)));
using uint64x4_t = uint64_t __attribute__((ext_vector_type(4)));
#else
using uint8x16_t [[gnu::vector_size(16)]] = uint8_t;
using int16x8_t [[gnu::vector_size(16)]] = int16_t;
using int32x4_t [[gnu::vector_size(16)]] = int32_t;
using int64x2_t [[gnu::vector_size(16)]] = int64_t;

using int16x16_t [[gnu::vector_size(32)]] = int16_t;
using int32x8_t [[gnu::vector_size(32)]] = int32_t;
using int64x4_t [[gnu::vector_size(32)]] = int64_t;

using uint8x32_t [[gnu::vector_size(32)]] = uint8_t;
using uint16x16_t [[gnu::vector_size(32)]] = uint16_t;
using uint32x8_t [[gnu::vector_size(32)]] = uint32_t;
using uint64x4_t [[gnu::vector_size(32)]] = uint64_t;
#endif

TEST(builtin_simd, builtin_simd)
{
    EXPECT_TRUE((std::is_same_v<detail::builtin_simd<int16_t, 8>::type, int16x8_t>));
    EXPECT_TRUE((std::is_same_v<detail::builtin_simd<int32_t, 4>::type, int32x4_t>));
    EXPECT_TRUE((std::is_same_v<detail::builtin_simd<int64_t, 2>::type, int64x2_t>));

    EXPECT_TRUE((std::is_same_v<detail::builtin_simd<uint16_t, 16>::type, uint16x16_t>));
    EXPECT_TRUE((std::is_same_v<detail::builtin_simd<uint32_t, 8>::type, uint32x8_t>));
    EXPECT_TRUE((std::is_same_v<detail::builtin_simd<uint64_t, 4>::type, uint64x4_t>));
}

TEST(builtin_simd, is_builtin_simd)
{
    EXPECT_FALSE(detail::is_builtin_simd<short>::value);
    EXPECT_FALSE(detail::is_builtin_simd<int>::value);
    EXPECT_FALSE(detail::is_builtin_simd<int[15]>::value);
    EXPECT_FALSE(detail::is_builtin_simd<int*>::value);

    EXPECT_TRUE(detail::is_builtin_simd<int16x8_t>::value);
    EXPECT_TRUE(detail::is_builtin_simd<int32x4_t>::value);
    EXPECT_TRUE(detail::is_builtin_simd<int64x2_t>::value);

    EXPECT_TRUE(detail::is_builtin_simd<uint16x16_t>::value);
    EXPECT_TRUE(detail::is_builtin_simd<uint32x8_t>::value);
    EXPECT_TRUE(detail::is_builtin_simd<uint64x4_t>::value);
}

TEST(builtin_simd, simd_traits)
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

    EXPECT_TRUE((std::is_same_v<simd_traits<int16x8_t>::mask_type, int16x8_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int32x4_t>::mask_type, int32x4_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int64x2_t>::mask_type, int64x2_t>));

    EXPECT_TRUE((std::is_same_v<simd_traits<int16x8_t>::swizzle_type, uint8x16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int32x4_t>::swizzle_type, uint8x16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<int64x2_t>::swizzle_type, uint8x16_t>));

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

    EXPECT_TRUE((std::is_same_v<simd_traits<uint16x16_t>::mask_type, int16x16_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint32x8_t>::mask_type, int32x8_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint64x4_t>::mask_type, int64x4_t>));

    EXPECT_TRUE((std::is_same_v<simd_traits<uint16x16_t>::swizzle_type, uint8x32_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint32x8_t>::swizzle_type, uint8x32_t>));
    EXPECT_TRUE((std::is_same_v<simd_traits<uint64x4_t>::swizzle_type, uint8x32_t>));
}

TEST(builtin_simd, default_simd_max_length)
{
    constexpr auto default_simd_max_length_v = detail::default_simd_max_length<detail::builtin_simd>;
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

TEST(builtin_simd, Simd)
{
    EXPECT_FALSE(Simd<short>);
    EXPECT_FALSE(Simd<int>);
    EXPECT_FALSE(Simd<int[15]>);
    EXPECT_FALSE(Simd<int*>);

    EXPECT_TRUE(Simd<int16x8_t>);
    EXPECT_TRUE(Simd<int32x4_t>);
    EXPECT_TRUE(Simd<int64x2_t>);
    EXPECT_TRUE(Simd<uint16x16_t>);
    EXPECT_TRUE(Simd<uint32x8_t>);
    EXPECT_TRUE(Simd<uint64x4_t>);
}
