// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
