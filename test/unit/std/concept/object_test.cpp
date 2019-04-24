// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <random>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>

#include "auxiliary.hpp"

using namespace seqan3;

TEST(object_concepts, Destructible)
{
    EXPECT_TRUE((std::Destructible<type_a>));
    EXPECT_TRUE((!std::Destructible<type_d>));
}

TEST(object_concepts, Constructible)
{
    EXPECT_TRUE((std::Constructible<type_a>));
    EXPECT_TRUE((std::Constructible<type_c, type_a>));
    EXPECT_TRUE((!std::Constructible<type_c, int>));
}

TEST(object_concepts, DefaultConstructible)
{
    EXPECT_TRUE((std::DefaultConstructible<type_a>));
    EXPECT_TRUE((!std::DefaultConstructible<type_d>));
}

TEST(object_concepts, MoveConstructible)
{
    EXPECT_TRUE((std::MoveConstructible<type_b>));
    EXPECT_TRUE((!std::MoveConstructible<type_d>));
}

TEST(object_concepts, CopyConstructible)
{
    EXPECT_TRUE((std::CopyConstructible<type_a>));
    EXPECT_TRUE((!std::CopyConstructible<type_b>));
}

TEST(object_concepts, Movable)
{
    EXPECT_TRUE((std::Movable<type_b>));
    EXPECT_TRUE((!std::Movable<type_d>));
}

TEST(object_concepts, Copyable)
{
    EXPECT_TRUE((std::Copyable<type_a>));
    EXPECT_TRUE((!std::Copyable<type_b>));
}

TEST(object_concepts, Semiregular)
{
    EXPECT_TRUE((std::Semiregular<type_a>));
    EXPECT_TRUE((std::Semiregular<type_c>));
    EXPECT_TRUE((!std::Semiregular<type_b>));
    EXPECT_TRUE((!std::Semiregular<type_d>));
}

TEST(object_concepts, Regular)
{
    EXPECT_TRUE((!std::Regular<type_a>));
    EXPECT_TRUE((!std::Regular<type_b>));
    EXPECT_TRUE((std::Regular<type_c>));
    EXPECT_TRUE((!std::Regular<type_d>));
}
