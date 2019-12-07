// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <forward_list>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/interleave.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/std/ranges>

#include <seqan3/test/pretty_printing.hpp>

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
    // explicitly call views::type_reduce
    auto v0 = views::type_reduce(u) | views::interleave(s, views::type_reduce(i));
    EXPECT_TRUE(ranges::equal(cmp, v0));
    EXPECT_EQ(cmpsize, v0.size());
    // don't call views::type_reduce
    auto v1 = u | views::interleave(s, i);
    EXPECT_TRUE(ranges::equal(cmp, v1));

    // function notation
    // explicitly call views::type_reduce
    auto v2{views::interleave(views::type_reduce(u), s, views::type_reduce(i))};
    EXPECT_TRUE(ranges::equal(cmp, v2));
    // don't call views::type_reduce
    auto v3{views::interleave(u, s, i)};
    EXPECT_TRUE(ranges::equal(cmp, v3));

    //combinability
    // explicitly call views::type_reduce
    auto v4 = views::type_reduce(u) | views::interleave(s, views::type_reduce(i)) | std::views::reverse | views::take(5);
    EXPECT_TRUE(ranges::equal(cmp_rev, v4));
    // don't call views::type_reduce
    auto v5 = u | views::interleave(s, i) | std::views::reverse | views::take(5);
    EXPECT_TRUE(ranges::equal(cmp_rev, v5));
}

TEST(view_interleave, concepts)
{
    // random_access_range, viewable_range, sized_range
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;
    auto v1 = detail::view_interleave(views::type_reduce(u), s, views::type_reduce(i));

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));

    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v1)>);

    // forward_range, viewable_range
    std::forward_list<dna4> u2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    dna4_vector i2{'G'_dna4};
    auto v2 = views::interleave(u2, s, i2);
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);

    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), dna4>));
}

#if 0 // blocked by https://github.com/ericniebler/range-v3/issues/1320
TEST(view_interleave, chunk_join)
{
    std::forward_list<dna4> u{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    dna4_vector i{'G'_dna4};
    size_t s = 2;

    dna4_vector cmp{"AAGAAGAA"_dna4};

    auto v1 = views::interleave(u, s, i);

    auto it = v1.begin();
    for (auto c : cmp)
    {
        EXPECT_EQ(c, *it);
        it++;
    }
}
#endif
