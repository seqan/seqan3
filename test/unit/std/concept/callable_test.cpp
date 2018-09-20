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

#include <gtest/gtest.h>

#include <random>

<<<<<<< HEAD
#include <seqan3/std/concept/callable.hpp>

#include "auxiliary.hpp"

using namespace seqan3;

TEST(invocable_concept, basic)
{
    EXPECT_TRUE((!invocable_concept<type_a, int, double,
                                    type_b>));
    EXPECT_TRUE((invocable_concept<std::random_device>));
    EXPECT_TRUE((invocable_concept<type_c, int, double,
                                   type_b>));
}

TEST(regular_invocable_concept, basic)
{
    EXPECT_TRUE((!regular_invocable_concept<type_a, int, double,
                                            type_b>));
//TODO(rrahn): Should not meet the regular_invocable_concept
//    EXPECT_TRUE((!regular_invocable_concept<std::random_device>));
    EXPECT_TRUE((regular_invocable_concept<type_c, int, double,
                                           type_b>));
}

TEST(predicate_concept, basic)
{
    EXPECT_TRUE((!predicate_concept<type_c, int, double,
                                    type_b>));
    EXPECT_TRUE((predicate_concept<type_b, int, double,
                                   type_b>));
}

TEST(relation_concept, basic)
{
    EXPECT_TRUE((!relation_concept<type_d, int, double>));
    EXPECT_TRUE((relation_concept<type_d, int, int>));
=======
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
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
}
