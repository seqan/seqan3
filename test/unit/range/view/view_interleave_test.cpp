// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <forward_list>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/interleave.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/ranges>

#include <gtest/gtest.h>

using namespace seqan3;

TEST(view_interleave, basic)
{
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;
    std::string cmp{"FOOinBARinBAXinBAT"};
    std::string cmp_rev{"TABni"};
    size_t cmpsize = 18;

    // pipe notation
    // explicitly call view::all
    auto v0 = view::all(u) | view::interleave(s, view::all(i));
    EXPECT_TRUE(ranges::equal(cmp, v0));
    EXPECT_EQ(cmpsize, v0.size());
    // don't call view::all
    auto v1 = u | view::interleave(s, i);
    EXPECT_TRUE(ranges::equal(cmp, v1));

    // function notation
    // explicitly call view::all
    auto v2{view::interleave(view::all(u), s, view::all(i))};
    EXPECT_TRUE(ranges::equal(cmp, v2));
    // don't call view::all
    auto v3{view::interleave(u, s, i)};
    EXPECT_TRUE(ranges::equal(cmp, v3));

    //combinability
    // explicitly call view::all
    auto v4 = view::all(u) | view::interleave(s, view::all(i)) | std::view::reverse | view::take(5);
    EXPECT_TRUE(ranges::equal(cmp_rev, v4));
    // don't call view::all
    auto v5 = u | view::interleave(s, i) | std::view::reverse | view::take(5);
    EXPECT_TRUE(ranges::equal(cmp_rev, v5));
}

TEST(view_interleave, concepts)
{
    // RandomAccessRange, ViewableRange, SizedRange
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;
    auto v1 = detail::view_interleave(view::all(u), s, view::all(i));

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), char>));

    EXPECT_FALSE(std::ranges::ContiguousRange<decltype(v1)>);

    // ForwardRange, ViewableRange
    std::forward_list<dna4> u2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    dna4_vector i2{'G'_dna4};
    auto v2 = view::interleave(u2, s, i2);
    EXPECT_TRUE(std::ranges::InputRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::View<decltype(v2)>);

    EXPECT_FALSE(std::ranges::ForwardRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::BidirectionalRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::RandomAccessRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::ContiguousRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::SizedRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v2)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v2), dna4>));
}

TEST(view_interleave, chunk_join)
{
    std::forward_list<dna4> u{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    dna4_vector i{'G'_dna4};
    size_t s = 2;

    dna4_vector cmp{"AAGAAGAA"_dna4};

    auto v1 = view::interleave(u, s, i);

    auto it = v1.begin();
    for (auto c : cmp)
    {
        EXPECT_EQ(c, *it);
        it++;
    }
}
