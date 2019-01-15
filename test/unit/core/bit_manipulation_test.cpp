// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/bit_manipulation.hpp>

using namespace seqan3::detail;

static constexpr size_t max_iterations = 1 << 15;

TEST(bit_manipulation, is_power_of_two)
{
    constexpr bool is_power_of_two0 = is_power_of_two(0);
    constexpr bool is_power_of_two1 = is_power_of_two(1);
    constexpr bool is_power_of_two2 = is_power_of_two(2);
    constexpr bool is_power_of_two3 = is_power_of_two(3);
    EXPECT_FALSE(is_power_of_two0);
    EXPECT_TRUE(is_power_of_two1);
    EXPECT_TRUE(is_power_of_two2);
    EXPECT_FALSE(is_power_of_two3);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_TRUE(is_power_of_two(power_of_two));

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_FALSE(is_power_of_two(i)) << i << " should not be a power of two.";
        }
    }
}

TEST(bit_manipulation, next_power_of_two)
{
    constexpr size_t next_power_of_two0 = next_power_of_two(0);
    constexpr size_t next_power_of_two1 = next_power_of_two(1);
    constexpr size_t next_power_of_two2 = next_power_of_two(2);
    constexpr size_t next_power_of_two3 = next_power_of_two(3);
    EXPECT_EQ(next_power_of_two0, 1u);
    EXPECT_EQ(next_power_of_two1, 1u);
    EXPECT_EQ(next_power_of_two2, 2u);
    EXPECT_EQ(next_power_of_two3, 4u);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_EQ(next_power_of_two(power_of_two), power_of_two);

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_EQ(next_power_of_two(i), next_power) << "The next power of two of " << i << " should be " << next_power;
        }
    }
}
