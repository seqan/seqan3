// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <iostream>
#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace seqan3::views
{
inline auto const deep_reverse = deep{std::views::reverse};
inline auto const deep_take = deep{std::views::take};
inline auto const deep_take2 = deep{std::views::take(2)};
} // namespace seqan3::views

using seqan3::operator""_dna5;

// ------------------------------------------------------------------
// no parameters
// ------------------------------------------------------------------

TEST(view_deep_reverse, basic)
{
    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    EXPECT_RANGE_EQ(foo | seqan3::views::deep{std::views::reverse}, "ATGCA"_dna5);

    // pipe notation
    EXPECT_RANGE_EQ(foo | seqan3::views::deep_reverse, "ATGCA"_dna5);

    // function notation
    EXPECT_RANGE_EQ(seqan3::views::deep_reverse(foo), "ATGCA"_dna5);

    // combinability
    EXPECT_RANGE_EQ(foo | seqan3::views::deep_reverse | std::views::reverse, "ACGTA"_dna5);
}

TEST(view_deep_reverse, deep)
{
    std::vector<seqan3::dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | seqan3::views::deep_reverse;

    ASSERT_EQ(std::ranges::size(v), 2u);
    EXPECT_RANGE_EQ(v[0], "ATGCA"_dna5);
    EXPECT_RANGE_EQ(v[1], "TACGT"_dna5);
}

TEST(view_deep_reverse, concepts)
{
    std::vector<seqan3::dna5_vector> vec{"ACGTA"_dna5, "TGCAT"_dna5};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::contiguous_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), seqan3::dna5_vector>));

    auto v1 = vec | seqan3::views::deep_reverse;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE(
        (std::ranges::output_range<decltype(v1), seqan3::dna5_vector>)); // view temporary returned in deep case

    auto v2 = v1 | std::views::reverse;
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v2)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_FALSE(
        (std::ranges::output_range<decltype(v2), seqan3::dna5_vector>)); // view temporary returned in deep case

    // v_elem has type std::ranges::reverse_view<std::ranges::ref_view<std::vector<seqan3::dna5> > >
    auto v_elem = v1[0];
    EXPECT_TRUE(std::ranges::input_range<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v_elem)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v_elem)>); // reverse_view drops contiguous_range
    EXPECT_TRUE(std::ranges::view<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v_elem)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v_elem)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v_elem), seqan3::dna5>));
}

// ------------------------------------------------------------------
// parameters preserved
// ------------------------------------------------------------------

TEST(view_deep_take, basic)
{
    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    EXPECT_RANGE_EQ(foo | seqan3::views::deep{std::views::take}(2), "AC"_dna5);

    // pipe notation
    EXPECT_RANGE_EQ(foo | seqan3::views::deep_take(2), "AC"_dna5);

    // function notation
    EXPECT_RANGE_EQ(seqan3::views::deep_take(foo, 2), "AC"_dna5);

    // combinability
    EXPECT_RANGE_EQ(foo | seqan3::views::deep_take(2) | std::views::reverse, "CA"_dna5);
}

TEST(view_deep_take, deep)
{
    std::vector<seqan3::dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5, "FOO"_dna5};

    // pipe notation
    auto v = foo | seqan3::views::deep_take(2);

    ASSERT_EQ(std::ranges::size(v), 3u);
    EXPECT_RANGE_EQ(v[0], "AC"_dna5);
    EXPECT_RANGE_EQ(v[1], "TG"_dna5);
    EXPECT_RANGE_EQ(v[2], "NN"_dna5);

    int i = 2;
    auto v2 = foo | seqan3::views::deep_take(i);

    ASSERT_EQ(std::ranges::size(v2), 3u);
    EXPECT_RANGE_EQ(v2[0], "AC"_dna5);
    EXPECT_RANGE_EQ(v2[1], "TG"_dna5);
    EXPECT_RANGE_EQ(v2[2], "NN"_dna5);

    // function notation
    auto v3 = seqan3::views::deep_take(foo, 2);

    ASSERT_EQ(std::ranges::size(v3), 3u);
    EXPECT_RANGE_EQ(v3[0], "AC"_dna5);
    EXPECT_RANGE_EQ(v3[1], "TG"_dna5);
    EXPECT_RANGE_EQ(v3[2], "NN"_dna5);
}

// ------------------------------------------------------------------
// parameters hardcoded
// ------------------------------------------------------------------

TEST(view_deep_take2, basic)
{
    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    EXPECT_RANGE_EQ(foo | seqan3::views::deep{std::views::take(2)}, "AC"_dna5);

    // pipe notation
    EXPECT_RANGE_EQ(foo | seqan3::views::deep_take2, "AC"_dna5);

    // function notation
    EXPECT_RANGE_EQ(seqan3::views::deep_take2(foo), "AC"_dna5);

    // combinability
    EXPECT_RANGE_EQ(foo | seqan3::views::deep_take2 | std::views::reverse, "CA"_dna5);
}

TEST(view_deep_take2, deep)
{
    std::vector<seqan3::dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5, "FOO"_dna5};

    // pipe notation
    auto v = foo | seqan3::views::deep_take2;

    ASSERT_EQ(std::ranges::size(v), 3u);
    EXPECT_RANGE_EQ(v[0], "AC"_dna5);
    EXPECT_RANGE_EQ(v[1], "TG"_dna5);
    EXPECT_RANGE_EQ(v[2], "NN"_dna5);

    // function notation
    auto v2 = seqan3::views::deep_take2(foo);

    ASSERT_EQ(std::ranges::size(v2), 3u);
    EXPECT_RANGE_EQ(v2[0], "AC"_dna5);
    EXPECT_RANGE_EQ(v2[1], "TG"_dna5);
    EXPECT_RANGE_EQ(v2[2], "NN"_dna5);
}
