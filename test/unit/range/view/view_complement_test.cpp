// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_complement, basic)
{
    dna5_vector foo{"ACGTA"_dna5};

    // pipe notation
    dna5_vector v = foo | view::complement;
    EXPECT_EQ(v, "TGCAT"_dna5);

    // function notation
    dna5_vector v2(view::complement(foo));
    EXPECT_EQ(v2, "TGCAT"_dna5);

    // combinability
    dna5_vector v3 = foo | view::complement | std::view::reverse;
    EXPECT_EQ(v3, "TACGT"_dna5);

    dna5_vector const bar{"ACGTA"_dna5};

    // const pipe notation
    dna5_vector v4 = bar | view::complement;
    EXPECT_EQ(v4, "TGCAT"_dna5);

    // const function notation
    dna5_vector v5(view::complement(bar));
    EXPECT_EQ(v5, "TGCAT"_dna5);

    // const combinability
    dna5_vector v6 = bar | view::complement | std::view::reverse;
    EXPECT_EQ(v6, "TACGT"_dna5);
}

TEST(view_complement, deep_view)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | view::complement;

    ASSERT_EQ(size(v), 2u);
    EXPECT_TRUE((std::ranges::equal(v[0], "TGCAT"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "ACGTA"_dna5)));

    std::vector<dna5_vector> const bar{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v2 = bar | view::complement;

    ASSERT_EQ(size(v2), 2u);
    EXPECT_TRUE((std::ranges::equal(v2[0], "TGCAT"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v2[1], "ACGTA"_dna5)));
}

TEST(view_complement, concepts)
{
    dna5_vector vec{"ACGTA"_dna5};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), dna5>));

    auto v1 = vec | view::complement;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), char>));

    dna5_vector const vec2{"ACGTA"_dna5};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec2)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec2)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec2)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(vec2), dna5>));

    auto v2 = vec2 | view::complement;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::View<decltype(v2)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v2)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v2)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v2), dna5>));
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v2), char>));
}
