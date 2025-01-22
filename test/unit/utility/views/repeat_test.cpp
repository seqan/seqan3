// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/views/repeat.hpp>

TEST(repeat_view, deduction_guide)
{
    int value = 0;
    int const & value_cref = value;

    seqan3::detail::repeat_view repeat_view1{value};
    EXPECT_SAME_TYPE(decltype(repeat_view1), seqan3::detail::repeat_view<int>);

    seqan3::detail::repeat_view repeat_view2{value_cref};
    EXPECT_SAME_TYPE(decltype(repeat_view2), seqan3::detail::repeat_view<int>);
}

TEST(general, construction)
{
    // char
    char chr{'A'};
    auto v = seqan3::views::repeat(chr);

    EXPECT_TRUE((std::is_default_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(v)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(v)>));

    // char const
    char const chr_const{'A'};
    auto v_const = seqan3::views::repeat(chr_const);

    EXPECT_TRUE((std::is_default_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(v_const)>));
}

TEST(general, concept)
{
    char chr{'A'};
    auto v = seqan3::views::repeat(chr);

    EXPECT_TRUE((std::ranges::range<decltype(v)>));
    EXPECT_TRUE((std::ranges::input_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::forward_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::bidirectional_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::random_access_range<decltype(v)>));
    EXPECT_FALSE((std::ranges::contiguous_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::view<decltype(v)>));
    EXPECT_FALSE((std::ranges::sized_range<decltype(v)>));
    EXPECT_FALSE((std::ranges::common_range<decltype(v)>));
    EXPECT_TRUE((std::ranges::output_range<decltype(v), char>));
}

TEST(general, iterator)
{
    auto v = seqan3::views::repeat('A');

    EXPECT_TRUE(v.begin() == v.begin());

    EXPECT_FALSE(v.begin() == v.end());
    EXPECT_FALSE(v.end() == v.begin());
    EXPECT_TRUE(v.begin() != v.end());
    EXPECT_TRUE(v.end() != v.begin());

    EXPECT_FALSE(std::ranges::cbegin(v) == std::ranges::cend(v));
    EXPECT_FALSE(std::ranges::cend(v) == std::ranges::cbegin(v));
    EXPECT_TRUE(std::ranges::cbegin(v) != std::ranges::cend(v));
    EXPECT_TRUE(std::ranges::cend(v) != std::ranges::cbegin(v));

    auto it = v.begin();
    EXPECT_EQ(*it, 'A');

    // random access iterator
    ++it;
    EXPECT_EQ(*it, 'A');
    it++;
    EXPECT_EQ(*it, 'A');
    --it;
    EXPECT_EQ(*it, 'A');
    it--;
    EXPECT_EQ(*it, 'A');
    it = it + 1;
    EXPECT_EQ(*it, 'A');
    it = 1 + it;
    EXPECT_EQ(*it, 'A');
    it += 1;
    EXPECT_EQ(*it, 'A');
    it = it - 1;
    EXPECT_EQ(*it, 'A');
    it = 1 - it;
    EXPECT_EQ(*it, 'A');
    it -= 1;
    EXPECT_EQ(*it, 'A');
    auto diff = (v.begin() + 1) - v.begin();
    EXPECT_EQ(diff, 1);

    // assignment
    *it = 'X';
    EXPECT_EQ(*it, 'X');

    auto cit = std::ranges::cbegin(v);
    cit = v.begin(); // assign it from const it
}

TEST(general, subscript_operator)
{
    auto v = seqan3::views::repeat('A');

    EXPECT_EQ(v[0], 'A');
    EXPECT_EQ(v[126], 'A');
    EXPECT_EQ(v[78'634'126], 'A');

    v[234] = 'X';

    EXPECT_EQ(v[0], 'X');
    EXPECT_EQ(v[126], 'X');
    EXPECT_EQ(v[78'634'126], 'X');
}

TEST(view, factory)
{
    // const char
    {
        char const chr{'X'};
        auto v = seqan3::views::repeat(chr);
        EXPECT_EQ(*v.begin(), chr);
    }

    // string
    {
        std::string str{"foobar"};
        auto v = seqan3::views::repeat(str);
        EXPECT_EQ(*v.begin(), str);
        EXPECT_EQ(v[2345], str);
    }

    // view
    {
        std::string const input{"foobar"};
        auto view = std::views::take(input, 3);
        auto v = seqan3::views::repeat(view);
        EXPECT_RANGE_EQ(*v.begin(), view);
    }

    // combinability
    {
        std::string str{"foobar"};
        auto v = seqan3::views::repeat(str) | std::views::take(3);

        using it_t = std::ranges::iterator_t<decltype(v)>;

        // Note: std::views::take doesn't handle infinite std::ranges::random_access_range's that are not
        //       std::ranges::sized_range, so we can't use std::ranges::size() to compute the size.
        // See https://eel.is/c++draft/range.take.view
        it_t it = v.begin();
        it_t sentinel = std::ranges::next(v.begin(), v.end()); // unfortunately linear in time
        auto size = sentinel - it;
        EXPECT_EQ(*it, str);
        EXPECT_EQ(size, 3u);
    }
}

constexpr char constexpr_class_and_iterator()
{
    auto v = seqan3::views::repeat('A');

    auto it = v.begin();
    ++it;
    it[234] = 'X';

    return *it;
}

constexpr char constexpr_view()
{
    char chr{'A'};
    auto v = seqan3::views::repeat(chr);
    v[1324] = 'X';

    return *v.begin();
}

TEST(general, constexpr_context)
{
    {
        constexpr char val = constexpr_class_and_iterator();
        EXPECT_EQ(val, 'X');
    }

    {
        constexpr char val = constexpr_view();
        EXPECT_EQ(val, 'X');
    }
}
