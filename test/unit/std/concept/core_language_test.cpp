// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
