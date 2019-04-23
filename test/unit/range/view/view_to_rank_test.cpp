// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_to_rank, basic)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    std::vector<unsigned> cmp{0,1,4,4,4,2,0,4,0};

    // pipe notation
    std::vector<unsigned> v = vec | view::to_rank;
    EXPECT_EQ(cmp, v);

    // function notation
    std::vector<unsigned> v2(view::to_rank(vec));
    EXPECT_EQ(cmp, v2);

    // combinability
    std::vector<unsigned> cmp2{0, 4, 0, 2, 4, 4, 4, 1, 0};
    std::vector<unsigned> v3 = vec | view::to_rank | std::view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_to_rank, concepts)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(ConstIterableRange<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), dna5>));

    auto v1 = vec | view::to_rank;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(ConstIterableRange<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), unsigned>));
}
