// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>

using seqan3::operator""_dna5;

TEST(view_char_strictly_to, basic)
{
    std::string vec{"ACTTTGATA"};
    seqan3::dna5_vector cmp{"ACTTTGATA"_dna5};

    // pipe notation
    EXPECT_RANGE_EQ(cmp, vec | seqan3::views::char_strictly_to<seqan3::dna5>);

    // function notation
    EXPECT_RANGE_EQ(cmp, seqan3::views::char_strictly_to<seqan3::dna5>(vec));

    // combinability
    seqan3::dna5_vector cmp2{"ATAGTTTCA"_dna5};
    EXPECT_RANGE_EQ(cmp2, vec | seqan3::views::char_strictly_to<seqan3::dna5> | std::views::reverse);
}

TEST(view_char_strictly_to, deep_view)
{
    std::vector<std::string> foo{"ACGTA", "TGCAT"};

    auto v = foo | seqan3::views::char_strictly_to<seqan3::dna5>;

    ASSERT_EQ(std::ranges::size(v), 2u);
    EXPECT_RANGE_EQ(v[0], "ACGTA"_dna5);
    EXPECT_RANGE_EQ(v[1], "TGCAT"_dna5);
}

TEST(view_char_strictly_to, concepts)
{
    std::string vec{"ACTTTGATA"};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::contiguous_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), char>));

    auto v1 = vec | seqan3::views::char_strictly_to<seqan3::dna5>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), seqan3::dna5>));
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));
}

TEST(view_char_strictly_to, exception)
{
    std::string foo = "ACGPTA";

    auto v = foo | seqan3::views::char_strictly_to<seqan3::dna5>;
    EXPECT_THROW(static_cast<void>(std::ranges::equal(v, "ACGNTA"_dna5)), seqan3::invalid_char_assignment);
}
