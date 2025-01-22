// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/views/convert.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

TEST(view_convert, basic)
{
    std::vector<int> vec{7, 5, 0, 5, 0, 0, 4, 8, -3};
    std::vector<bool> cmp{1, 1, 0, 1, 0, 0, 1, 1, 1};

    // pipe notation
    EXPECT_RANGE_EQ(cmp, vec | seqan3::views::convert<bool>);

    // function notation
    EXPECT_RANGE_EQ(cmp, seqan3::views::convert<bool>(vec));

    // combinability
    std::vector<bool> cmp2{1, 1, 1, 0, 0, 1, 0, 1, 1};
    EXPECT_RANGE_EQ(cmp2, vec | seqan3::views::convert<bool> | std::views::reverse);
}

TEST(view_convert, explicit_conversion)
{
    seqan3::dna5_vector vec{"ACGNTNGGN"_dna5};
    seqan3::dna4_vector cmp{"ACGATAGGA"_dna4};

    // pipe notation
    EXPECT_RANGE_EQ(cmp, vec | seqan3::views::convert<seqan3::dna4>);

    // function notation
    EXPECT_RANGE_EQ(cmp, seqan3::views::convert<seqan3::dna4>(vec));

    // combinability
    seqan3::dna4_vector cmp2{"AGGATAGCA"_dna4};
    EXPECT_RANGE_EQ(cmp2, vec | seqan3::views::convert<seqan3::dna4> | std::views::reverse);
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
