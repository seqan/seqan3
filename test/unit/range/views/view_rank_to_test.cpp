// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/rank_to.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_rank_to, basic)
{
    std::vector<unsigned> vec{0,1,4,4,4,2,0,4,0};
    dna5_vector cmp{"ACTTTGATA"_dna5};

    // pipe notation
    dna5_vector v = vec | views::rank_to<dna5> | views::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    dna5_vector v2(views::rank_to<dna5>(vec) | views::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    dna5_vector cmp2{"ATAGTTTCA"_dna5};
    dna5_vector v3 = vec | views::rank_to<dna5> | std::views::reverse | views::to<std::vector>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_rank_to, concepts)
{
    std::vector<unsigned> vec{0,1,3,3,3,2,0,3,0};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), unsigned>));

    auto v1 = vec | views::rank_to<dna5>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), unsigned>));
}
