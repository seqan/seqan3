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
#include <seqan3/std/concept/comparison.hpp>
=======
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

#include "auxiliary.hpp"

using namespace seqan3;

<<<<<<< HEAD
TEST(weakly_equality_comparable_with_concept, basic)
{
    EXPECT_TRUE((weakly_equality_comparable_with_concept<type_a, type_b>));
    EXPECT_TRUE((!weakly_equality_comparable_with_concept<type_a, type_c>));
}

TEST(equality_comparable_with_concept, basic)
{
    EXPECT_TRUE((!equality_comparable_with_concept<type_a, type_b>));
    EXPECT_TRUE((equality_comparable_with_concept<type_b, type_d>));
}

TEST(strict_totally_ordered_with_concept, basic)
{
    EXPECT_TRUE((!strict_totally_ordered_with_concept<type_a, type_b>));
    EXPECT_TRUE((strict_totally_ordered_with_concept<type_b, type_d>));
=======
TEST(comparison_concepts, WeaklyEqualityComparableWith)
{
    EXPECT_TRUE((std::detail::WeaklyEqualityComparableWith<type_a, type_b>));
    EXPECT_TRUE((!std::detail::WeaklyEqualityComparableWith<type_a, type_c>));
}

TEST(comparison_concepts, EqualityComparableWith)
{
    EXPECT_TRUE((!std::EqualityComparableWith<type_a, type_b>));
    EXPECT_TRUE((std::EqualityComparableWith<type_b, type_d>));
}

TEST(comparison_concepts, StrictTotallyOrderedWith)
{
    EXPECT_TRUE((!std::StrictTotallyOrderedWith<type_a, type_b>));
    EXPECT_TRUE((std::StrictTotallyOrderedWith<type_b, type_d>));
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
}
