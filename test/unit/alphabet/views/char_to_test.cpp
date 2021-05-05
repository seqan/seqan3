// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <iostream>
#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>

using seqan3::operator""_dna5;

TEST(view_char_to, basic)
{
    std::string vec{"ACTTTGATA"};
    seqan3::dna5_vector cmp{"ACTTTGATA"_dna5};

    // pipe notation
    EXPECT_RANGE_EQ(cmp, vec | seqan3::views::char_to<seqan3::dna5>);

    // function notation
    EXPECT_RANGE_EQ(cmp, seqan3::views::char_to<seqan3::dna5>(vec));

    // combinability
    seqan3::dna5_vector cmp2{"ATAGTTTCA"_dna5};
    EXPECT_RANGE_EQ(cmp2, vec | seqan3::views::char_to<seqan3::dna5> | std::views::reverse);
}

TEST(view_char_to, deep_view)
{
    std::vector<std::string> foo{"ACGTA", "TGCAT"};

    auto v = foo | seqan3::views::char_to<seqan3::dna5>;

    ASSERT_EQ(std::ranges::size(v), 2u);
    EXPECT_RANGE_EQ(v[0], "ACGTA"_dna5);
    EXPECT_RANGE_EQ(v[1], "TGCAT"_dna5);
}

TEST(view_char_to, concepts)
{
    std::string vec{"ACTTTGATA"};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), char>));

    auto v1 = vec | seqan3::views::char_to<seqan3::dna5>;
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
}
