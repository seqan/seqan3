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

TEST(comparison_concepts, weakly_equality_comparable_with)
{
    EXPECT_TRUE((seqan3::detail::weakly_equality_comparable_with<type_a, type_b>));
    EXPECT_TRUE((!seqan3::detail::weakly_equality_comparable_with<type_a, type_c>));
}

TEST(comparison_concepts, equality_comparable)
{
    EXPECT_TRUE((!std::equality_comparable_with<type_a, type_b>));
    EXPECT_TRUE((std::equality_comparable_with<type_b, type_d>));
}

TEST(comparison_concepts, totally_ordered)
{
    EXPECT_TRUE((!std::totally_ordered_with<type_a, type_b>));
    EXPECT_TRUE((std::totally_ordered_with<type_b, type_d>));
}
