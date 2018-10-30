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
    EXPECT_EQ(next_power_of_two0, 0);
    EXPECT_EQ(next_power_of_two1, 1);
    EXPECT_EQ(next_power_of_two2, 2);
    EXPECT_EQ(next_power_of_two3, 4);

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
