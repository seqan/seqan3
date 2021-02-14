// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/concepts>
#include <random>

#include <seqan3/utility/detail/exposition_only_concept.hpp>

#include "auxiliary.hpp"

TEST(object_concepts, destructible)
{
    EXPECT_TRUE((std::destructible<type_a>));
    EXPECT_TRUE((!std::destructible<type_d>));
}

TEST(object_concepts, constructible_from)
{
    EXPECT_TRUE((std::constructible_from<type_a>));
    EXPECT_TRUE((std::constructible_from<type_c, type_a>));
    EXPECT_TRUE((!std::constructible_from<type_c, int>));
}

TEST(object_concepts, default_initializable)
{
    EXPECT_TRUE((std::default_initializable<type_a>));
    EXPECT_TRUE((!std::default_initializable<type_d>));
}

TEST(object_concepts, move_constructible)
{
    EXPECT_TRUE((std::move_constructible<type_b>));
    EXPECT_TRUE((!std::move_constructible<type_d>));
}

TEST(object_concepts, copy_constructible)
{
    EXPECT_TRUE((std::copy_constructible<type_a>));
    EXPECT_TRUE((!std::copy_constructible<type_b>));
}

TEST(object_concepts, movable)
{
    EXPECT_TRUE((std::movable<type_b>));
    EXPECT_TRUE((!std::movable<type_d>));
}

TEST(object_concepts, copyable)
{
    EXPECT_TRUE((std::copyable<type_a>));
    EXPECT_TRUE((!std::copyable<type_b>));
}

TEST(object_concepts, semiregular)
{
    EXPECT_TRUE((std::semiregular<type_a>));
    EXPECT_TRUE((std::semiregular<type_c>));
    EXPECT_TRUE((!std::semiregular<type_b>));
    EXPECT_TRUE((!std::semiregular<type_d>));
}

TEST(object_concepts, regular)
{
    EXPECT_TRUE((!std::regular<type_a>));
    EXPECT_TRUE((!std::regular<type_b>));
    EXPECT_TRUE((std::regular<type_c>));
    EXPECT_TRUE((!std::regular<type_d>));
}
