// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <random>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>

#include "auxiliary.hpp"

using namespace seqan3;

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
}
