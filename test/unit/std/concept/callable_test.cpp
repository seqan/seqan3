// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <random>

#include <seqan3/std/concepts>

#include "auxiliary.hpp"

TEST(callable_concepts, invocable)
{
    EXPECT_TRUE((!std::invocable<type_a, int, double, type_b>));
    EXPECT_TRUE((std::invocable<std::random_device>));
    EXPECT_TRUE((std::invocable<type_c, int, double, type_b>));
}

TEST(callable_concepts, regular_invocable)
{
    EXPECT_TRUE((!std::regular_invocable<type_a, int, double, type_b>));
//TODO(rrahn): Should not meet the std::regular_invocable
//    EXPECT_TRUE((!std::regular_invocable<std::random_device>));
    EXPECT_TRUE((std::regular_invocable<type_c, int, double, type_b>));
}

TEST(callable_concepts, predicate)
{
    EXPECT_TRUE((!std::predicate<type_c, int, double, type_b>));
    EXPECT_TRUE((std::predicate<type_b, int, double, type_b>));
}

TEST(callable_concepts, relation)
{
    EXPECT_TRUE((!std::relation<type_d, int, double>));
    EXPECT_TRUE((std::relation<type_d, int, int>));
}
