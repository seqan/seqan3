// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gap.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

INSTANTIATE_TYPED_TEST_CASE_P(gap, alphabet, gap);
INSTANTIATE_TYPED_TEST_CASE_P(gap, alphabet_constexpr, gap);

TEST(gap_test, default_initialization)
{
    gap gap1;
    gap gap2{};
    gap gap3 = gap{};

    EXPECT_EQ(gap1.to_rank(), 0);
    EXPECT_EQ(gap2.to_rank(), 0);
    EXPECT_EQ(gap3.to_rank(), 0);
    EXPECT_EQ(gap1.to_char(), '-');
    EXPECT_EQ(gap2.to_char(), '-');
    EXPECT_EQ(gap3.to_char(), '-');
}

TEST(gap_test, relations)
{
    EXPECT_EQ(gap{}, gap{});
    EXPECT_LE(gap{}, gap{});
    EXPECT_GE(gap{}, gap{});
}

TEST(gap_test, assign_char)
{
    EXPECT_EQ(gap{}.assign_char('-'), gap{});
    EXPECT_EQ(gap{}.assign_char('x'), gap{});
}

TEST(gap_test, to_rank)
{
    EXPECT_EQ(gap{}.to_rank(), 0);
}

TEST(gap_test, assign_rank)
{
    EXPECT_EQ(gap{}.assign_rank(0), gap{});
    // EXPECT_EQ(gap{}.assign_rank(13), gap{});
}
