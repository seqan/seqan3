// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <type_traits>

#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/detail/builtin_simd.hpp>
#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

#ifdef __clang__
using int8x16_t = int8_t __attribute__((ext_vector_type(16)));
using int16x8_t = int16_t __attribute__((ext_vector_type(8)));
using int32x4_t = int32_t __attribute__((ext_vector_type(4)));
using int64x2_t = int64_t __attribute__((ext_vector_type(2)));

using uint8x16_t = uint8_t __attribute__((ext_vector_type(16)));
using uint16x8_t = uint16_t __attribute__((ext_vector_type(8)));
using uint32x4_t = uint32_t __attribute__((ext_vector_type(4)));
using uint64x2_t = uint64_t __attribute__((ext_vector_type(2)));

using int8x32_t = int8_t __attribute__((ext_vector_type(32)));
using int16x16_t = int16_t __attribute__((ext_vector_type(16)));
using int32x8_t = int32_t __attribute__((ext_vector_type(8)));
using int64x4_t = int64_t __attribute__((ext_vector_type(4)));

using uint8x32_t = uint8_t __attribute__((ext_vector_type(32)));
using uint16x16_t = uint16_t __attribute__((ext_vector_type(16)));
using uint32x8_t = uint32_t __attribute__((ext_vector_type(8)));
using uint64x4_t = uint64_t __attribute__((ext_vector_type(4)));
#else
using int8x16_t [[gnu::vector_size(16)]] = int8_t;
using int16x8_t [[gnu::vector_size(16)]] = int16_t;
using int32x4_t [[gnu::vector_size(16)]] = int32_t;
using int64x2_t [[gnu::vector_size(16)]] = int64_t;

using uint8x16_t [[gnu::vector_size(16)]] = uint8_t;
using uint16x8_t [[gnu::vector_size(16)]] = uint16_t;
using uint32x4_t [[gnu::vector_size(16)]] = uint32_t;
using uint64x2_t [[gnu::vector_size(16)]] = uint64_t;

using int8x32_t [[gnu::vector_size(32)]] = int8_t;
using int16x16_t [[gnu::vector_size(32)]] = int16_t;
using int32x8_t [[gnu::vector_size(32)]] = int32_t;
using int64x4_t [[gnu::vector_size(32)]] = int64_t;

using uint8x32_t [[gnu::vector_size(32)]] = uint8_t;
using uint16x16_t [[gnu::vector_size(32)]] = uint16_t;
using uint32x8_t [[gnu::vector_size(32)]] = uint32_t;
using uint64x4_t [[gnu::vector_size(32)]] = uint64_t;
#endif

// Wrapping `type` and `template_type` in a namespace triggered different bugs
namespace incomplete
{

struct type;

template <typename t>
struct template_type;

} // namespace incomplete

// Types that have a subscript operator
template <typename type>
using subscript_types = seqan3::type_list<type[15],
                                          type const[15],
                                          type[15][15],
                                          type const[15][15],
                                          type *,
                                          type const *,
                                          type * [15][15],
                                          type **,
                                          type const **,
                                          type const * const *,
                                          type ** [15][15],
                                          type const ** [15][15],
                                          type const * const * [15][15],
                                          type ***,
                                          type const ***,
                                          type const * const **,
                                          type const * const * const *,
                                          type *** [15][15],
                                          type const *** [15][15],
                                          type const * const ** [15][15],
                                          type const * const * const * [15][15],
                                          type const * const * const * const[15][15]>;

TEST(builtin_simd, builtin_simd)
{
    EXPECT_TRUE((std::is_same_v<seqan3::detail::builtin_simd<int16_t, 8>::type, int16x8_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::builtin_simd<int32_t, 4>::type, int32x4_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::builtin_simd<int64_t, 2>::type, int64x2_t>));

    EXPECT_TRUE((std::is_same_v<seqan3::detail::builtin_simd<uint16_t, 16>::type, uint16x16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::builtin_simd<uint32_t, 8>::type, uint32x8_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::builtin_simd<uint64_t, 4>::type, uint64x4_t>));
}

TEST(builtin_simd, is_builtin_simd)
{
    EXPECT_FALSE(seqan3::detail::is_builtin_simd<short>::value);
    EXPECT_FALSE(seqan3::detail::is_builtin_simd<int>::value);
    EXPECT_FALSE(seqan3::detail::is_builtin_simd<incomplete::type>::value);
    EXPECT_FALSE(seqan3::detail::is_builtin_simd<incomplete::template_type<int>>::value);

    auto is_not_builtin_simd = [](auto id)
    {
        using type = typename decltype(id)::type;
        EXPECT_FALSE(seqan3::detail::is_builtin_simd<type>::value);
    };

    seqan3::detail::for_each<subscript_types<short>>(is_not_builtin_simd);
    seqan3::detail::for_each<subscript_types<int>>(is_not_builtin_simd);
    seqan3::detail::for_each<subscript_types<incomplete::type>>(is_not_builtin_simd);
    seqan3::detail::for_each<subscript_types<incomplete::template_type<int>>>(is_not_builtin_simd);

    EXPECT_TRUE(seqan3::detail::is_builtin_simd<int16x8_t>::value);
    EXPECT_TRUE(seqan3::detail::is_builtin_simd<int32x4_t>::value);
    EXPECT_TRUE(seqan3::detail::is_builtin_simd<int64x2_t>::value);

    EXPECT_TRUE(seqan3::detail::is_builtin_simd<uint16x16_t>::value);
    EXPECT_TRUE(seqan3::detail::is_builtin_simd<uint32x8_t>::value);
    EXPECT_TRUE(seqan3::detail::is_builtin_simd<uint64x4_t>::value);
}

TEST(builtin_simd, simd_traits)
{
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int16x8_t>::scalar_type, int16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x4_t>::scalar_type, int32_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int64x2_t>::scalar_type, int64_t>));

    EXPECT_EQ(seqan3::simd::simd_traits<int16x8_t>::length, 8u);
    EXPECT_EQ(seqan3::simd::simd_traits<int32x4_t>::length, 4u);
    EXPECT_EQ(seqan3::simd::simd_traits<int64x2_t>::length, 2u);

    EXPECT_EQ(seqan3::simd::simd_traits<int16x8_t>::max_length, 16u);
    EXPECT_EQ(seqan3::simd::simd_traits<int32x4_t>::max_length, 16u);
    EXPECT_EQ(seqan3::simd::simd_traits<int64x2_t>::max_length, 16u);

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int16x8_t>::mask_type,
                                decltype(std::declval<int16x8_t>() == std::declval<int16x8_t>())>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x4_t>::mask_type,
                                decltype(std::declval<int32x4_t>() == std::declval<int32x4_t>())>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int64x2_t>::mask_type,
                                decltype(std::declval<int64x2_t>() == std::declval<int64x2_t>())>));

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int16x8_t>::swizzle_type, uint8x16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x4_t>::swizzle_type, uint8x16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int64x2_t>::swizzle_type, uint8x16_t>));

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int8x16_t>::rebind<uint8_t>, uint8x16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int16x8_t>::rebind<uint16_t>, uint16x8_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x4_t>::rebind<uint32_t>, uint32x4_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int64x2_t>::rebind<uint64_t>, uint64x2_t>));

    // avx2 (256bit)

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint16x16_t>::scalar_type, uint16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint32x8_t>::scalar_type, uint32_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint64x4_t>::scalar_type, uint64_t>));

    EXPECT_EQ(seqan3::simd::simd_traits<uint16x16_t>::length, 16u);
    EXPECT_EQ(seqan3::simd::simd_traits<uint32x8_t>::length, 8u);
    EXPECT_EQ(seqan3::simd::simd_traits<uint64x4_t>::length, 4u);

    EXPECT_EQ(seqan3::simd::simd_traits<uint16x16_t>::max_length, 32u);
    EXPECT_EQ(seqan3::simd::simd_traits<uint32x8_t>::max_length, 32u);
    EXPECT_EQ(seqan3::simd::simd_traits<uint64x4_t>::max_length, 32u);

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint16x16_t>::mask_type,
                                decltype(std::declval<int16x16_t>() == std::declval<int16x16_t>())>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint32x8_t>::mask_type,
                                decltype(std::declval<uint32x8_t>() == std::declval<uint32x8_t>())>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint64x4_t>::mask_type,
                                decltype(std::declval<uint64x4_t>() == std::declval<uint64x4_t>())>));

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint16x16_t>::swizzle_type, uint8x32_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint32x8_t>::swizzle_type, uint8x32_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<uint64x4_t>::swizzle_type, uint8x32_t>));

    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int8x32_t>::rebind<uint8_t>, uint8x32_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int16x16_t>::rebind<uint16_t>, uint16x16_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x8_t>::rebind<uint32_t>, uint32x8_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int64x4_t>::rebind<uint64_t>, uint64x4_t>));

    // Currently the simd_concept does not support floating point scalar types.
    // EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<float>::rebind<uint32_t>, uint32x4_t>));
    // EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<float>::rebind<uint32_t>, uint32x8_t>));
    // EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x4_t>::rebind<float>, float>));
    // EXPECT_TRUE((std::is_same_v<seqan3::simd::simd_traits<int32x8_t>::rebind<float>, float>));
}

TEST(builtin_simd, default_simd_max_length)
{
    constexpr auto default_simd_max_length_v = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;
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

TEST(builtin_simd, simd)
{
    EXPECT_FALSE(seqan3::simd::simd_concept<short>);
    EXPECT_FALSE(seqan3::simd::simd_concept<int>);
    EXPECT_FALSE(seqan3::simd::simd_concept<incomplete::type>);
    EXPECT_FALSE(seqan3::simd::simd_concept<incomplete::template_type<int>>);

    auto fails_simd = [](auto id)
    {
        using type = typename decltype(id)::type;
        EXPECT_FALSE(seqan3::simd::simd_concept<type>);
    };

    seqan3::detail::for_each<subscript_types<short>>(fails_simd);
    seqan3::detail::for_each<subscript_types<int>>(fails_simd);
    seqan3::detail::for_each<subscript_types<incomplete::type>>(fails_simd);
    seqan3::detail::for_each<subscript_types<incomplete::template_type<int>>>(fails_simd);

    EXPECT_TRUE(seqan3::simd::simd_concept<int16x8_t>);
    EXPECT_TRUE(seqan3::simd::simd_concept<int32x4_t>);
    EXPECT_TRUE(seqan3::simd::simd_concept<int64x2_t>);
    EXPECT_TRUE(seqan3::simd::simd_concept<uint16x16_t>);
    EXPECT_TRUE(seqan3::simd::simd_concept<uint32x8_t>);
    EXPECT_TRUE(seqan3::simd::simd_concept<uint64x4_t>);
}
