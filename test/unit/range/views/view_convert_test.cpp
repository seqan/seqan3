// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/ranges>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

TEST(view_convert, basic)
{
    std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};
    std::vector<bool> cmp{1, 1, 0, 1, 0, 0, 1, 1, 1};

    // pipe notation
    std::vector<bool> v = vec | seqan3::views::convert<bool> | seqan3::views::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    std::vector<bool> v2(seqan3::views::convert<bool>(vec) | seqan3::views::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    std::vector<bool> cmp2{1, 1, 1, 0, 0, 1, 0, 1, 1};
    std::vector<bool> v3 = vec | seqan3::views::convert<bool> | std::views::reverse | seqan3::views::to<std::vector>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, explicit_conversion)
{
    seqan3::dna5_vector vec{"ACGNTNGGN"_dna5};
    seqan3::dna4_vector cmp{"ACGATAGGA"_dna4};

    // pipe notation
    seqan3::dna4_vector v = vec | seqan3::views::convert<seqan3::dna4> | seqan3::views::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    seqan3::dna4_vector v2(seqan3::views::convert<seqan3::dna4>(vec) | seqan3::views::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    seqan3::dna4_vector cmp2{"AGGATAGCA"_dna4};
    seqan3::dna4_vector v3 = vec
                           | seqan3::views::convert<seqan3::dna4>
                           | std::views::reverse
                           | seqan3::views::to<std::vector>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, concepts)
{
    seqan3::dna5_vector vec{"ACGNTNGGN"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), seqan3::dna5>));

    auto v1 = vec | seqan3::views::convert<seqan3::dna4>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), seqan3::dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), seqan3::dna4>));
}
