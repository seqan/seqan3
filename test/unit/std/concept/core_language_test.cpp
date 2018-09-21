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

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/iterator>

#include "auxiliary.hpp"

using namespace seqan3;

TEST(core_language_concepts, Same)
{
    EXPECT_TRUE((std::Same<int, int>));
    EXPECT_TRUE((!std::Same<int, char>));
}

TEST(core_language_concepts, DerivedFrom)
{
    EXPECT_TRUE((std::DerivedFrom<type_b, type_a>));
    EXPECT_TRUE((!std::DerivedFrom<type_a, type_b>));
}

TEST(implicitly_convertible_to_concept, basic)
{
    EXPECT_TRUE((implicitly_convertible_to_concept<type_b, type_c>));
    EXPECT_TRUE((!implicitly_convertible_to_concept<type_c, type_b>));
    EXPECT_TRUE((!implicitly_convertible_to_concept<type_a, type_c>));
}

TEST(explicitly_convertible_to_concept, basic)
{
    EXPECT_TRUE((explicitly_convertible_to_concept<type_b, type_c>));
    EXPECT_TRUE((!explicitly_convertible_to_concept<type_c, type_b>));
    EXPECT_TRUE((explicitly_convertible_to_concept<type_a, type_c>));
}

TEST(core_language_concepts, ConvertibleTo)
{
    EXPECT_TRUE((std::ConvertibleTo<type_b, type_c>));
    EXPECT_TRUE((!std::ConvertibleTo<type_c, type_b>));
    EXPECT_TRUE((!std::ConvertibleTo<type_a, type_c>));
}

TEST(core_language_concepts, CommonReference)
{
    EXPECT_TRUE((std::CommonReference<int32_t, int16_t>));
    EXPECT_TRUE((!std::CommonReference<int32_t, type_c>));
}

TEST(core_language_concepts, Common)
{
    EXPECT_TRUE((std::Common<type_a, type_b>));
    EXPECT_TRUE((!std::Common<type_a, type_c>));
}

TEST(core_language_concepts, Integral)
{
    EXPECT_TRUE((std::Integral<int>));
    EXPECT_TRUE((!std::Integral<float>));
}

TEST(core_language_concepts, SignedIntegral)
{
    EXPECT_TRUE((std::SignedIntegral<int>));
    EXPECT_TRUE((!std::SignedIntegral<unsigned>));
}

TEST(core_language_concepts, UnsignedIntegral)
{
    EXPECT_TRUE((!std::UnsignedIntegral<int>));
    EXPECT_TRUE((std::UnsignedIntegral<unsigned>));
}

TEST(core_language_concepts, Assignable)
{
    EXPECT_TRUE((std::Assignable<type_a &, type_a const &>));
    EXPECT_TRUE((std::Assignable<type_c &, type_b const &>));
    EXPECT_TRUE((!std::Assignable<type_a &, type_c &>));
}

TEST(core_language_concepts, Swappable)
{
    EXPECT_TRUE((std::Swappable<type_a>));
    EXPECT_TRUE((std::Swappable<type_b>));
}

TEST(core_language_concepts, SwappableWith)
{
    EXPECT_TRUE((std::SwappableWith<type_a &, type_a &>));
    EXPECT_TRUE((!std::SwappableWith<type_b, type_c>));
}
