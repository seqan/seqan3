// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using int16x_t = seqan3::simd::simd_type<int16_t>::type;

TEST(simd, auto_length)
{
    if constexpr(simd_traits<int16x_t>::max_length == 64u)
    {
        using int16x32_t = seqan3::simd::simd_type<int16_t, 32>::type; // avx512 512bit
        EXPECT_TRUE((std::is_same_v<int16x_t, int16x32_t>));
        EXPECT_EQ(simd_traits<int16x_t>::length, 32u);
    }
    else if constexpr(simd_traits<int16x_t>::max_length == 32u)
    {
        using int16x16_t = seqan3::simd::simd_type<int16_t, 16>::type; // avx2 256bit
        EXPECT_TRUE((std::is_same_v<int16x_t, int16x16_t>));
        EXPECT_EQ(simd_traits<int16x_t>::length, 16u);
    }
    else if constexpr(simd_traits<int16x_t>::max_length == 16u)
    {
        using int16x8_t = seqan3::simd::simd_type<int16_t, 8>::type; // sse4 128bit
        EXPECT_TRUE((std::is_same_v<int16x_t, int16x8_t>));
        EXPECT_EQ(simd_traits<int16x_t>::length, 8u);
    }
    else if constexpr(simd_traits<int16x_t>::max_length == 1u)
    {
        using int16x1_t = seqan3::simd::simd_type<int16_t, 1>::type; // scalar fallback
        EXPECT_TRUE((std::is_same_v<int16x_t, int16x1_t>));
        EXPECT_EQ(simd_traits<int16x_t>::length, 1u);
    }
    else
    {
        ADD_FAILURE() << "unknown max_length: " << simd_traits<int16x_t>::max_length;
    }
}

TEST(simd, standard_construction)
{
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<int16x_t>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<int16x_t>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<int16x_t>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<int16x_t>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<int16x_t>));
    EXPECT_TRUE((std::is_nothrow_destructible_v<int16x_t>));
    EXPECT_TRUE((std::is_nothrow_swappable_v<int16x_t>));
}

template <simd::simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 1>)
{
    simd = simd_t{a};
}

template <simd::simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 8>)
{
    simd = simd_t{a, a, a, a, a, a, a, a};
}

template <simd::simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 16>)
{
    simd = simd_t{a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a};
}

template <simd::simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 32>)
{
    simd = simd_t{a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a};
}

template <simd::simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 64>)
{
    simd = simd_t{a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a};
}

TEST(simd, construct)
{
    int16x_t a;
    construct_test(a, 4, std::integral_constant<size_t, simd_traits<int16x_t>::length>{});

    for (unsigned i = 0u; i < simd_traits<int16x_t>::length; ++i)
        EXPECT_EQ(a[i], 4);
}
