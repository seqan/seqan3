// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/as_const.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/to_lower.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

using seqan3::operator""_dna5;

TEST(view_as_const, basic)
{
    std::string vec{"ACTTTGATA"};

    // pipe notation
    auto v = vec | seqan3::views::as_const;
    EXPECT_TRUE((std::ranges::equal(vec, v)));

    // function notation
    auto v2(seqan3::views::as_const(vec));
    EXPECT_TRUE((std::ranges::equal(vec, v2)));

    // combinability
    seqan3::dna5_vector vec2{"ACGTA"_dna5};
    seqan3::dna5_vector v3 = vec2
                           | seqan3::views::complement
                           | seqan3::views::as_const
                           | seqan3::views::to<std::vector>;
    EXPECT_EQ("TGCAT"_dna5, v3);
}

TEST(view_as_const, concepts)
{
    std::string vec{"ACTTTGATA"};
    auto v1 = vec | seqan3::views::as_const;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));

    EXPECT_TRUE((std::is_same_v<decltype(v1[0]), char const &>));

    auto v2 = vec | seqan3::views::to_lower | seqan3::views::as_const; // to_lower generates values
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));

    EXPECT_TRUE((std::is_same_v<decltype(v2[0]), char>)); // don't add const-ness to values
}
