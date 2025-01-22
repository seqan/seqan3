// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <deque>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>
#include <seqan3/utility/views/slice.hpp>

// ============================================================================
//  view_slice
// ============================================================================

TEST(view_slice, regular)
{
    using namespace std::literals;

    std::string const vec{"foobar"};

    // pipe notation
    EXPECT_RANGE_EQ("oob"sv, vec | seqan3::views::slice(1, 4));

    // function notation
    EXPECT_RANGE_EQ("oob"sv, seqan3::views::slice(vec, 1, 4));

    // combinability
    EXPECT_RANGE_EQ("o"sv, vec | seqan3::views::slice(0, 4) | seqan3::views::slice(1, 3) | std::views::take(1));
    EXPECT_RANGE_EQ("abo"sv, vec | std::views::reverse | seqan3::views::slice(1, 4) | std::views::take(3));

    // store arg
    auto a0 = seqan3::views::slice(1, 4);
    EXPECT_RANGE_EQ("oob"sv, vec | a0);

    // store combined
    auto a1 = seqan3::views::slice(0, 4) | seqan3::views::slice(1, 3) | std::views::take(1);
    EXPECT_RANGE_EQ("o"sv, vec | a1);

    // store combined in middle
    auto a2 = std::views::reverse | seqan3::views::slice(1, 4) | std::views::take(3);
    EXPECT_RANGE_EQ("abo"sv, vec | a2);
}

TEST(view_slice, concepts)
{
    std::vector vec{1, 2, 3};

    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), int>));

    auto v1 = vec | seqan3::views::slice(1, 4);

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), int>));

    auto v2 = vec | seqan3::views::single_pass_input | seqan3::views::slice(1, 4);

    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(v2)>, false);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), int>)); // single_pass_input loses it
}

TEST(view_slice, underlying_is_shorter)
{
    using namespace std::literals;

    std::string vec{"foobar"};
    EXPECT_NO_THROW((seqan3::views::slice(vec, 1, 4))); // no parsing

    // full parsing on conversion
    EXPECT_RANGE_EQ("oob"sv, vec | seqan3::views::single_pass_input | seqan3::views::slice(1, 4));
}

TEST(view_slice, end_before_begin)
{
    std::string vec{"foobar"};
    EXPECT_THROW(seqan3::views::slice(vec, 4, 1), std::invalid_argument);
}

TEST(view_slice, type_erasure)
{
    { // string overload
        std::string const urange{"foobar"};

        auto v = seqan3::views::slice(urange, 1, 4);

        EXPECT_SAME_TYPE(decltype(v), std::string_view);
        EXPECT_RANGE_EQ(v, urange.substr(1, 3));
    }

    { // stringview overload
        std::string_view urange{"foobar"};

        auto v = seqan3::views::slice(urange, 1, 4);

        EXPECT_SAME_TYPE(decltype(v), std::string_view);
        EXPECT_RANGE_EQ(v, urange.substr(1, 3));
    }

    { // contiguous overload
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::slice(urange, 1, 4);

        EXPECT_SAME_TYPE(decltype(v), (std::span<int, std::dynamic_extent>));
        EXPECT_RANGE_EQ(v, (std::vector{2, 3, 4}));
    }

    { // contiguous overload
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::slice(urange, 1, 4);

        EXPECT_SAME_TYPE(decltype(v), (std::span<int, std::dynamic_extent>));
        EXPECT_RANGE_EQ(v, (std::vector{2, 3, 4}));
    }

    { // random-access overload
        std::deque<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::slice(urange, 1, 4);

        EXPECT_TRUE((std::same_as<
                     decltype(v),
                     std::ranges::subrange<typename std::deque<int>::iterator, typename std::deque<int>::iterator>>));
        EXPECT_RANGE_EQ(v, (std::vector{2, 3, 4}));
    }
}
