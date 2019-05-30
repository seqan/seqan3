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
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_to_char, basic)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    std::string cmp{"ACTTTGATA"};

    // pipe notation
    std::string v = vec | view::to_char | std::ranges::to<std::string>;
    EXPECT_EQ(cmp, v);

    // function notation
    std::string v2(view::to_char(vec) | std::ranges::to<std::string>);
    EXPECT_EQ(cmp, v2);

    // combinability
    std::string cmp2{"ATAGTTTCA"};
    std::string v3 = vec | view::to_char | std::view::reverse | std::ranges::to<std::string>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_to_char, concepts)
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

    auto v1 = vec | view::to_char;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(ConstIterableRange<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), char>));
}
