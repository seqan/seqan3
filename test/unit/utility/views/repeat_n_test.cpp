// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/repeat_n.hpp>

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
        std::string str{"foobar"};
        auto view = str | std::views::take(3);
        auto v = seqan3::views::repeat_n(view, 5);
        EXPECT_RANGE_EQ(*v.begin(), std::string{"foo"});
    }

    // combinability
    {
        std::string str{"foobar"};
        auto v = seqan3::views::repeat_n(str, 2)
               | std::views::transform(
                     [](auto & str)
                     {
                         return str.substr(3);
                     });
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
