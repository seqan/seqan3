// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
<<<<<<< HEAD
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
=======
// Copyright (C) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (C) 2016-2018, Knut Reinert & MPI Molekulare Genetik
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
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

#include <gtest/gtest.h>

#include <random>

<<<<<<< HEAD
#include <seqan3/std/concept/object.hpp>
=======
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

#include "auxiliary.hpp"

using namespace seqan3;

<<<<<<< HEAD
TEST(destructible_concept, basic)
{
    EXPECT_TRUE((destructible_concept<type_a>));
    EXPECT_TRUE((!destructible_concept<type_d>));
}

TEST(constructible_concept, basic)
{
    EXPECT_TRUE((constructible_concept<type_a>));
    EXPECT_TRUE((constructible_concept<type_c,
                                       type_a>));
    EXPECT_TRUE((!constructible_concept<type_c,
                                        type_a,
                                        type_b>));
}

TEST(default_constructible_concept, basic)
{
    EXPECT_TRUE((default_constructible_concept<type_a>));
    EXPECT_TRUE((!default_constructible_concept<type_d>));
}

TEST(move_constructible_concept, basic)
{
    EXPECT_TRUE((move_constructible_concept<type_b>));
    EXPECT_TRUE((!move_constructible_concept<type_d>));
}

TEST(copy_constructible_concept, basic)
{
    EXPECT_TRUE((copy_constructible_concept<type_a>));
    EXPECT_TRUE((!copy_constructible_concept<type_b>));
}

TEST(movable_concept, basic)
{
    EXPECT_TRUE((movable_concept<type_b>));
    EXPECT_TRUE((!movable_concept<type_d>));
}

TEST(copyable_concept, basic)
{
    EXPECT_TRUE((copyable_concept<type_a>));
    EXPECT_TRUE((!copyable_concept<type_b>));
}

TEST(semi_regular_concept, basic)
{
    EXPECT_TRUE((semi_regular_concept<type_a>));
    EXPECT_TRUE((semi_regular_concept<type_c>));
    EXPECT_TRUE((!semi_regular_concept<type_b>));
    EXPECT_TRUE((!semi_regular_concept<type_d>));
}

TEST(regular_concept, basic)
{
    EXPECT_TRUE((!regular_concept<type_a>));
    EXPECT_TRUE((!regular_concept<type_b>));
    EXPECT_TRUE((regular_concept<type_c>));
    EXPECT_TRUE((!regular_concept<type_d>));
=======
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
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
}
