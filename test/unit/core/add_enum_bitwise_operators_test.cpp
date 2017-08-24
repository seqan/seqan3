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

#include <seqan3/core/add_enum_bitwise_operators.hpp>

using namespace seqan3;

enum class my_enum
{
    ZERO = 0,
    VAL1 = 1,
    VAL2 = 2,
    COMB = 3
};

template <>
constexpr bool seqan3::add_enum_bitwise_operators<my_enum> = true;

TEST(add_enum_bitwise_operators, AND)
{
    my_enum e = my_enum::VAL1;
    my_enum e2 = e & my_enum::VAL2;
    EXPECT_EQ(e2, my_enum::ZERO);
}

TEST(add_enum_bitwise_operators, OR)
{
    my_enum e = my_enum::VAL1;
    my_enum e2 = e | my_enum::VAL2;
    EXPECT_EQ(e2, my_enum::COMB);
}

TEST(add_enum_bitwise_operators, XOR)
{
    my_enum e = my_enum::VAL1;
    my_enum e2 = e ^ my_enum::VAL2;
    EXPECT_EQ(e2, my_enum::COMB);
}

TEST(add_enum_bitwise_operators, NOT)
{
    my_enum e = my_enum::VAL1;
    my_enum e2 = ~e;
    EXPECT_NE(e, e2);
    e2 = ~e2;
    EXPECT_EQ(e, e2);
}

TEST(add_enum_bitwise_operators, AND_ASSIGN)
{
    my_enum e = my_enum::VAL1;
    e &= my_enum::VAL2;
    EXPECT_EQ(e, my_enum::ZERO);
}

TEST(add_enum_bitwise_operators, OR_ASSIGN)
{
    my_enum e = my_enum::VAL1;
    e |= my_enum::VAL2;
    EXPECT_EQ(e, my_enum::COMB);
}

TEST(add_enum_bitwise_operators, XOR_ASSIGN)
{
    my_enum e = my_enum::VAL1;
    e ^= my_enum::VAL2;
    EXPECT_EQ(e, my_enum::COMB);
}
