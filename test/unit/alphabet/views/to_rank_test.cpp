// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/to_rank.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>

using seqan3::operator""_dna5;

TEST(view_to_rank, basic)
{
    seqan3::dna5_vector vec{"ACTTTGATA"_dna5};
    std::vector<uint8_t> cmp{0,1,4,4,4,2,0,4,0};

    // pipe notation
    EXPECT_RANGE_EQ(cmp, vec | seqan3::views::to_rank);

    // function notation
    EXPECT_RANGE_EQ(cmp, seqan3::views::to_rank(vec));

    // combinability
    std::vector<uint8_t> cmp2{0, 4, 0, 2, 4, 4, 4, 1, 0};
    EXPECT_RANGE_EQ(cmp2, vec | seqan3::views::to_rank | std::views::reverse);
}

TEST(view_to_rank, concepts)
{
    seqan3::dna5_vector vec{"ACTTTGATA"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), seqan3::dna5>));

    auto v1 = vec | seqan3::views::to_rank;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), seqan3::dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), uint8_t>));
}
