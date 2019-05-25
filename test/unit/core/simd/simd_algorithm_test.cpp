// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/test/simd_utility.hpp>

#include <iostream>
template <typename t>
class simd_algorithm_generic : public ::testing::Test
{
public:
};

using test_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t>;

TYPED_TEST_CASE(simd_algorithm_generic, test_types);

using namespace seqan3;

TEST(simd_algorithm, fill)
{
    using simd_type = simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < simd_traits<simd_type>::length; ++i)
        expect[i] = 4;

    constexpr simd_type result = fill<simd_type>(4);
    SIMD_EQ(result, expect);
}

TEST(simd_algorithm, iota)
{
    using simd_type = simd_type_t<int16_t, 8>;

    simd_type expect{};
    for (size_t i = 0; i < simd_traits<simd_type>::length; ++i)
        expect[i] = i;

    constexpr simd_type result = iota<simd_type>(0);
    SIMD_EQ(result, expect);
}

TYPED_TEST(simd_algorithm_generic, unpack_hi)
{
    using simd_t = typename simd_type<TypeParam>::type;

    // Nothing to test if the simd type is just scalar.
    // TODO: report gcc9.1 fails for intrinsic _mm512_permutex2var_epi8 on mac: g++-9 (Homebrew GCC 9.1.0) 9.1.0
    if constexpr (simd_traits<simd_t>::length == 1 || simd_traits<simd_t>::length == 64)
    {
        return;
    }
    else
    {
        simd_t lhs = iota<simd_t>(1);
        simd_t rhs = iota<simd_t>(simd_traits<simd_t>::length + 1);

        simd_t res = unpack_hi(lhs, rhs);

        simd_t cmp{};

        uint32_t j = simd_traits<simd_t>::length / 2;
        for (uint32_t i = 0; i < simd_traits<simd_t>::length; i += 2, ++j)
        {
            cmp[i]     = lhs[j];
            cmp[i + 1] = rhs[j];
        }

        SIMD_EQ(res, cmp);
    }
}

TYPED_TEST(simd_algorithm_generic, unpack_lo)
{
    using simd_t = typename simd_type<TypeParam>::type;

    // Nothing to test if the simd type is just scalar.
    // TODO: report gcc9.1 fails for intrinsic _mm512_permutex2var_epi8 on mac: g++-9 (Homebrew GCC 9.1.0) 9.1.0
    if constexpr (simd_traits<simd_t>::length == 1 || simd_traits<simd_t>::length == 64)
    {
        return;
    }
    else
    {

        simd_t lhs = iota<simd_t>(1);
        simd_t rhs = iota<simd_t>(simd_traits<simd_t>::length + 1);

        simd_t res = unpack_lo(lhs, rhs);

        simd_t cmp{};

        for (uint32_t i = 0, j = 0; i < simd_traits<simd_t>::length; i += 2, ++j)
        {
            cmp[i]     = lhs[j];
            cmp[i + 1] = rhs[j];
        }

        SIMD_EQ(res, cmp);
    }
}
