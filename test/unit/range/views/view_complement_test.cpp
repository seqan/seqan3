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
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

using seqan3::operator""_dna5;

TEST(view_complement, basic)
{
    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation
    seqan3::dna5_vector v = foo | seqan3::views::complement | seqan3::views::to<std::vector>;
    EXPECT_EQ(v, "TGCAT"_dna5);

    // function notation
    seqan3::dna5_vector v2(seqan3::views::complement(foo) | seqan3::views::to<std::vector>);
    EXPECT_EQ(v2, "TGCAT"_dna5);

    // combinability
    seqan3::dna5_vector v3 = foo | seqan3::views::complement | std::views::reverse | seqan3::views::to<std::vector>;
    EXPECT_EQ(v3, "TACGT"_dna5);

    seqan3::dna5_vector const bar{"ACGTA"_dna5};

    // const pipe notation
    seqan3::dna5_vector v4 = bar | seqan3::views::complement | seqan3::views::to<std::vector>;
    EXPECT_EQ(v4, "TGCAT"_dna5);

    // const function notation
    seqan3::dna5_vector v5(seqan3::views::complement(bar) | seqan3::views::to<std::vector>);
    EXPECT_EQ(v5, "TGCAT"_dna5);

    // const combinability
    seqan3::dna5_vector v6 = bar | seqan3::views::complement | std::views::reverse | seqan3::views::to<std::vector>;
    EXPECT_EQ(v6, "TACGT"_dna5);
}

TEST(view_complement, deep_view)
{
    std::vector<seqan3::dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | seqan3::views::complement;

    ASSERT_EQ(size(v), 2u);
    EXPECT_TRUE((std::ranges::equal(v[0], "TGCAT"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "ACGTA"_dna5)));

    std::vector<seqan3::dna5_vector> const bar{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v2 = bar | seqan3::views::complement;

    ASSERT_EQ(size(v2), 2u);
    EXPECT_TRUE((std::ranges::equal(v2[0], "TGCAT"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v2[1], "ACGTA"_dna5)));
}

TEST(view_complement, concepts)
{
    seqan3::dna5_vector vec{"ACGTA"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), seqan3::dna5>));

    auto v1 = vec | seqan3::views::complement;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), seqan3::dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));

    seqan3::dna5_vector const vec2{"ACGTA"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec2)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec2)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec2)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(vec2), seqan3::dna5>));

    auto v2 = vec2 | seqan3::views::complement;
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v2)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), seqan3::dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), char>));
}
