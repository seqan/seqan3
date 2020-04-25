// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/type_traits/iterator.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/test/expect_range_eq.hpp>

TEST(general, construction)
{
    // char
    char chr{'A'};
    auto v = seqan3::views::repeat_n(chr, 4);

    EXPECT_TRUE((std::is_default_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(v)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(v)>));

    // char const
    char const chr_const{'A'};
    auto v_const = seqan3::views::repeat_n(chr_const, 20);

    EXPECT_TRUE((std::is_default_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(v_const)>));
}

TEST(general, concept)
{
    char chr{'A'};
    auto v = seqan3::views::repeat_n(chr, 10);

    EXPECT_TRUE((std::ranges::range<decltype(v)>));
    EXPECT_TRUE((std::ranges::input_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::forward_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::bidirectional_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::random_access_range<decltype(v)>));
    EXPECT_FALSE((std::ranges::contiguous_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::view<decltype(v)>));
    EXPECT_TRUE((std::ranges::sized_range<decltype(v)>));
    EXPECT_FALSE((std::ranges::common_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::output_range<decltype(v), char>));
}

TEST(view, factory)
{
    // const char
    {
        char const chr{'X'};
        auto v = seqan3::views::repeat_n(chr, 3);
        EXPECT_EQ(v.size(), 3u);
        EXPECT_RANGE_EQ(v, (std::vector<char>{chr, chr, chr}));
    }

    // string
    {
        std::string str{"foobar"};
        auto v = seqan3::views::repeat_n(str, 2);
        EXPECT_EQ(v.size(), 2u);
        EXPECT_EQ(*v.begin(), str);
        EXPECT_EQ(v[0], str);
    }

    // view
    {
        auto view = std::string{"foobar"} | seqan3::views::persist | std::views::take(3);
        auto v = seqan3::views::repeat_n(view, 5);
        EXPECT_RANGE_EQ(*v.begin(), std::string{"foo"});
    }

    // combinability
    {
        std::string str{"foobar"};
        auto v = seqan3::views::repeat_n(str, 2) | std::views::transform([] (auto & str) { return str.substr(3); });
        EXPECT_RANGE_EQ(v, (std::vector<std::string>{"bar", "bar"}));
    }
}

constexpr char constexpr_view()
{
    char chr{'A'};
    auto v = seqan3::views::repeat_n(chr, 10);
    v[0] = 'X';

    return *v.begin();
}

TEST(general, constexpr_context)
{
    constexpr char val = constexpr_view();
    EXPECT_EQ(val, 'X');
}
