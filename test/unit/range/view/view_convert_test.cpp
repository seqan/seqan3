// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_convert, basic)
{
    std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};
    std::vector<bool> cmp{1, 1, 0, 1, 0, 0, 1, 1, 1};

    // pipe notation
    std::vector<bool> v = vec | view::convert<bool> | std::ranges::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    std::vector<bool> v2(view::convert<bool>(vec) | std::ranges::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    std::vector<bool> cmp2{1, 1, 1, 0, 0, 1, 0, 1, 1};
    std::vector<bool> v3 = vec | view::convert<bool> | std::view::reverse | std::ranges::to<std::vector>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, explicit_conversion)
{
    dna5_vector vec{"ACGNTNGGN"_dna5};
    dna4_vector cmp{"ACGATAGGA"_dna4};

    // pipe notation
    dna4_vector v = vec | view::convert<dna4> | std::ranges::to<std::vector>;
    EXPECT_EQ(cmp, v);

    // function notation
    dna4_vector v2(view::convert<dna4>(vec) | std::ranges::to<std::vector>);
    EXPECT_EQ(cmp, v2);

    // combinability
    dna4_vector cmp2{"AGGATAGCA"_dna4};
    dna4_vector v3 = vec | view::convert<dna4> | std::view::reverse | std::ranges::to<std::vector>;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, concepts)
{
    dna5_vector vec{"ACGNTNGGN"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), dna5>));

    auto v1 = vec | view::convert<dna4>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), dna4>));
}
