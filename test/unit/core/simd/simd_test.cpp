#include <gtest/gtest.h>

#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

#if defined(__SSE4_2__) || defined(__AVX2__) || defined(__AVX512F__)

using simd_type = seqan3::simd<int16_t>::type;

TEST(simd, auto_length)
{
    if constexpr(simd_traits<simd_type>::max_length == 64)
    {
        using int16x32_t = seqan3::simd<int16_t, 32>::type; // avx512 512bit
        EXPECT_TRUE((std::is_same_v<simd_type, int16x32_t>));
        EXPECT_EQ(simd_traits<simd_type>::length, 32);
    }
    else if constexpr(simd_traits<simd_type>::max_length == 32)
    {
        using int16x16_t = seqan3::simd<int16_t, 16>::type; // avx2 256bit
        EXPECT_TRUE((std::is_same_v<simd_type, int16x16_t>));
        EXPECT_EQ(simd_traits<simd_type>::length, 16);
    }
    else if constexpr(simd_traits<simd_type>::max_length == 16)
    {
        using int16x8_t = seqan3::simd<int16_t, 8>::type; // sse4 128bit
        EXPECT_TRUE((std::is_same_v<simd_type, int16x8_t>));
        EXPECT_EQ(simd_traits<simd_type>::length, 8);
    }
}

TEST(simd, standard_construction)
{
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<simd_type>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<simd_type>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<simd_type>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<simd_type>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<simd_type>));
    EXPECT_TRUE((std::is_nothrow_destructible_v<simd_type>));
    EXPECT_TRUE((std::is_nothrow_swappable_v<simd_type>));
}

template <simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 8>)
{
    simd = simd_t{a, a, a, a, a, a, a, a};
}

template <simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 16>)
{
    simd = simd_t{a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a};
}

template <simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 32>)
{
    simd = simd_t{a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a};
}

template <simd_concept simd_t>
void construct_test(simd_t & simd, typename simd_traits<simd_t>::scalar_type a, std::integral_constant<size_t, 64>)
{
    simd = simd_t{a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a,
                  a, a, a, a, a, a, a, a, a, a, a, a, a, a, a, a};
}

TEST(simd, construct)
{
    simd_type a;
    construct_test(a, 4, std::integral_constant<size_t, simd_traits<simd_type>::length>{});

    for (unsigned i = 0u; i < simd_traits<simd_type>::length; ++i)
        EXPECT_EQ(a[i], 4);
}

#endif
