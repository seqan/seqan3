// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <random>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/iterator>

#include "auxiliary.hpp"

TEST(core_language_concepts, same_as)
{
    EXPECT_TRUE((std::same_as<int, int>));
    EXPECT_TRUE((!std::same_as<int, char>));
}

TEST(core_language_concepts, derived_from)
{
    EXPECT_TRUE((std::derived_from<type_b, type_a>));
    EXPECT_TRUE((!std::derived_from<type_a, type_b>));
}

TEST(implicitly_convertible_to, basic)
{
    EXPECT_TRUE((seqan3::implicitly_convertible_to<type_b, type_c>));
    EXPECT_TRUE((!seqan3::implicitly_convertible_to<type_c, type_b>));
    EXPECT_TRUE((!seqan3::implicitly_convertible_to<type_a, type_c>));
}

TEST(explicitly_convertible_to, basic)
{
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_b, type_c>));
    EXPECT_TRUE((!seqan3::explicitly_convertible_to<type_c, type_b>));
    EXPECT_TRUE((seqan3::explicitly_convertible_to<type_a, type_c>));
}

TEST(core_language_concepts, convertible_to)
{
    EXPECT_TRUE((std::convertible_to<type_b, type_c>));
    EXPECT_TRUE((!std::convertible_to<type_c, type_b>));
    EXPECT_TRUE((!std::convertible_to<type_a, type_c>));
}

TEST(core_language_concepts, common_reference_with)
{
    EXPECT_TRUE((std::common_reference_with<int32_t, int16_t>));
    EXPECT_TRUE((!std::common_reference_with<int32_t, type_c>));
}

TEST(core_language_concepts, common_with)
{
    EXPECT_TRUE((std::common_with<type_a, type_b>));
    EXPECT_TRUE((!std::common_with<type_a, type_c>));
}

TEST(core_language_concepts, integral)
{
    EXPECT_TRUE((std::integral<int>));
    EXPECT_TRUE((!std::integral<float>));
}

TEST(core_language_concepts, signed_integral)
{
    EXPECT_TRUE((std::signed_integral<int>));
    EXPECT_TRUE((!std::signed_integral<unsigned>));
}

TEST(core_language_concepts, unsigned_integral)
{
    EXPECT_TRUE((!std::unsigned_integral<int>));
    EXPECT_TRUE((std::unsigned_integral<unsigned>));
}

TEST(core_language_concepts, assignable_from)
{
    EXPECT_TRUE((std::assignable_from<type_a &, type_a const &>));
    EXPECT_TRUE((std::assignable_from<type_c &, type_b const &>));
    EXPECT_TRUE((!std::assignable_from<type_a &, type_c &>));
}

TEST(core_language_concepts, swappable)
{
    EXPECT_TRUE((std::swappable<type_a>));
    EXPECT_TRUE((std::swappable<type_b>));
}

TEST(core_language_concepts, swappable_with)
{
    EXPECT_TRUE((std::swappable_with<type_a &, type_a &>));
    EXPECT_TRUE((!std::swappable_with<type_b, type_c>));
}
