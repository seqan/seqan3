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

TEST(view_rank_to, basic)
{
    using seqan3::operator""_dna5;

    std::vector<unsigned> vec{0,1,4,4,4,2,0,4,0};
    seqan3::dna5_vector cmp{"ACTTTGATA"_dna5};

    // pipe notation
    seqan3::dna5_vector v = vec | seqan3::views::rank_to<seqan3::dna5> | seqan3::views::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    seqan3::dna5_vector v2(seqan3::views::rank_to<seqan3::dna5>(vec) | seqan3::views::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    seqan3::dna5_vector cmp2{"ATAGTTTCA"_dna5};
    seqan3::dna5_vector v3 = vec
                           | seqan3::views::rank_to<seqan3::dna5>
                           | std::views::reverse
                           | seqan3::views::to<std::vector>;
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
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), unsigned>));

    auto v1 = vec | seqan3::views::rank_to<seqan3::dna5>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), seqan3::dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), unsigned>));
}
