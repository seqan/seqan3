// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/test/simd_utility.hpp>

#include <algorithm>
#include <numeric>
#include <iostream>

using namespace seqan3;

template <typename simd_t>
simd_t transform_iota(int offset)
{
    return detail::simd_transform<simd_t>([offset] (size_t const i)
    {
        return offset + i;
    });
}

template <simd_concept simd_t>
constexpr simd_t transform_iota_constexpr(int offset)
{
    return detail::simd_transform_constexpr<simd_t>([offset] (size_t const i)
    {
        return offset + i;
    });
}

TEST(simd_transform_constexpr, nullary_iota)
{
    using simd_type = simd_type_t<int16_t, 8>;
    constexpr size_t length = simd_traits<simd_type>::length;

    constexpr simd_type result = transform_iota_constexpr<simd_type>(4);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = 4 + i;

    SIMD_EQ(result, expect);
}

#ifndef __clang__
TEST(simd_transform_constexpr, unary_add)
{
    using simd_type = simd_type_t<int16_t, 8>;
    constexpr size_t length = simd_traits<simd_type>::length;

    constexpr simd_type four_iota = transform_iota_constexpr<simd_type>(4);
    constexpr simd_type result = detail::simd_transform_constexpr<simd_type>([] (size_t, auto four_iota_i)
    {
        return four_iota_i + 6;
    }, four_iota);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = i + 4 + 6;

    SIMD_EQ(result, expect);
}

TEST(simd_transform_constexpr, binary_multiply)
{
    using simd_type = simd_type_t<int16_t, 8>;
    constexpr size_t length = simd_traits<simd_type>::length;

    constexpr simd_type four_iota = transform_iota_constexpr<simd_type>(4);
    constexpr simd_type two_iota = transform_iota_constexpr<simd_type>(2);
    constexpr simd_type result = detail::simd_transform_constexpr<simd_type>([] (size_t, auto four_iota_i, auto two_iota_i)
    {
        return four_iota_i * two_iota_i;
    }, four_iota, two_iota);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = (i + 4) * (i + 2);

    SIMD_EQ(result, expect);
}

TEST(simd_transform_constexpr, ternary_blend)
{
    using simd_type = simd_type_t<int16_t, 8>;
    using mask_type = typename simd_traits<simd_type>::mask_type;
    constexpr size_t length = simd_traits<simd_type>::length;

    constexpr simd_type four_iota = transform_iota_constexpr<simd_type>(4);
    constexpr simd_type two_iota = transform_iota_constexpr<simd_type>(2);
    constexpr mask_type mask = detail::simd_transform_constexpr<mask_type>([] (size_t i)
    {
        return (i % 3);
    });

    constexpr simd_type result = detail::simd_transform_constexpr<simd_type>([] (size_t, auto four_iota_i, auto two_iota_i, auto mask_i)
    {
        return mask_i ? four_iota_i : two_iota_i;
    }, four_iota, two_iota, mask);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = (i % 3) ? (i + 4) : (i + 2);

    SIMD_EQ(result, expect);
}
#endif

TEST(simd_transform, nullary_iota)
{
    using simd_type = simd_type_t<int16_t, 8>;
    constexpr size_t length = simd_traits<simd_type>::length;

    simd_type result = transform_iota<simd_type>(4);
    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = 4 + i;

    SIMD_EQ(result, expect);
}

TEST(simd_transform, unary_add)
{
    using simd_type = simd_type_t<int16_t, 8>;
    constexpr size_t length = simd_traits<simd_type>::length;

    simd_type four_iota = transform_iota<simd_type>(4);
    simd_type result = detail::simd_transform<simd_type>([] (size_t, auto four_iota_i)
    {
        return four_iota_i + 6;
    }, four_iota);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = i + 4 + 6;

    SIMD_EQ(result, expect);
}

TEST(simd_transform, binary_multiply)
{
    using simd_type = simd_type_t<int16_t, 8>;
    constexpr size_t length = simd_traits<simd_type>::length;

    simd_type four_iota = transform_iota<simd_type>(4);
    simd_type two_iota = transform_iota<simd_type>(2);
    simd_type result = detail::simd_transform<simd_type>([] (size_t, auto four_iota_i, auto two_iota_i)
    {
        return four_iota_i * two_iota_i;
    }, four_iota, two_iota);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = (i + 4) * (i + 2);

    SIMD_EQ(result, expect);
}

TEST(simd_transform, ternary_blend)
{
    using simd_type = simd_type_t<int16_t, 8>;
    using mask_type = typename simd_traits<simd_type>::mask_type;
    constexpr size_t length = simd_traits<simd_type>::length;

    simd_type four_iota = transform_iota<simd_type>(4);
    simd_type two_iota = transform_iota<simd_type>(2);
    mask_type mask = detail::simd_transform<mask_type>([] (size_t i)
    {
        return (i % 3);
    });

    simd_type result = detail::simd_transform<simd_type>([] (size_t, auto four_iota_i, auto two_iota_i, auto mask_i)
    {
        return mask_i ? four_iota_i : two_iota_i;
    }, four_iota, two_iota, mask);

    simd_type expect{};
    for (size_t i = 0; i < length; ++i)
        expect[i] = (i % 3) ? (i + 4) : (i + 2);

    SIMD_EQ(result, expect);
}
