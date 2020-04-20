// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gap.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(gap, alphabet_, seqan3::gap, );
INSTANTIATE_TYPED_TEST_SUITE_P(gap, semi_alphabet_test, seqan3::gap, );
INSTANTIATE_TYPED_TEST_SUITE_P(gap, alphabet_constexpr, seqan3::gap, );
INSTANTIATE_TYPED_TEST_SUITE_P(gap, semi_alphabet_constexpr, seqan3::gap, );

TEST(gap_test, default_initialization)
{
    seqan3::gap gap1;
    seqan3::gap gap2{};
    seqan3::gap gap3 = seqan3::gap{};

    EXPECT_EQ(gap1.to_rank(), 0);
    EXPECT_EQ(gap2.to_rank(), 0);
    EXPECT_EQ(gap3.to_rank(), 0);
    EXPECT_EQ(gap1.to_char(), '-');
    EXPECT_EQ(gap2.to_char(), '-');
    EXPECT_EQ(gap3.to_char(), '-');
}

TEST(gap_test, relations)
{
    EXPECT_EQ(seqan3::gap{}, seqan3::gap{});
    EXPECT_LE(seqan3::gap{}, seqan3::gap{});
    EXPECT_GE(seqan3::gap{}, seqan3::gap{});
}

TEST(gap_test, assign_char)
{
    EXPECT_EQ(seqan3::gap{}.assign_char('-'), seqan3::gap{});
    EXPECT_EQ(seqan3::gap{}.assign_char('x'), seqan3::gap{});
}

TEST(gap_test, to_rank)
{
    EXPECT_EQ(seqan3::gap{}.to_rank(), 0);
}

TEST(gap_test, assign_rank)
{
    EXPECT_EQ(seqan3::gap{}.assign_rank(0), seqan3::gap{});
    // EXPECT_EQ(seqan3::gap{}.assign_rank(13), seqan3::gap{});
}
