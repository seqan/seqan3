// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <random>

#include <seqan3/std/concepts>

#include "auxiliary.hpp"

TEST(callable_concepts, Invocable)
{
    EXPECT_TRUE((!std::Invocable<type_a, int, double, type_b>));
    EXPECT_TRUE((std::Invocable<std::random_device>));
    EXPECT_TRUE((std::Invocable<type_c, int, double, type_b>));
}

TEST(callable_concepts, RegularInvocable)
{
    EXPECT_TRUE((!std::RegularInvocable<type_a, int, double, type_b>));
//TODO(rrahn): Should not meet the std::RegularInvocable
//    EXPECT_TRUE((!std::RegularInvocable<std::random_device>));
    EXPECT_TRUE((std::RegularInvocable<type_c, int, double, type_b>));
}

TEST(callable_concepts, Predicate)
{
    EXPECT_TRUE((!std::Predicate<type_c, int, double, type_b>));
    EXPECT_TRUE((std::Predicate<type_b, int, double, type_b>));
}

TEST(callable_concepts, Relation)
{
    EXPECT_TRUE((!std::Relation<type_d, int, double>));
    EXPECT_TRUE((std::Relation<type_d, int, int>));
}
