// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/type_traits/iterator.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/repeat.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(general, construction)
{
    // char
    char chr{'A'};
    auto v = view::repeat(chr);

    EXPECT_TRUE((std::is_default_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(v)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(v)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(v)>));

    // char const
    char const chr_const{'A'};
    auto v_const = view::repeat(chr_const);

    EXPECT_TRUE((std::is_default_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(v_const)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(v_const)>));
}

TEST(general, concept)
{
    char chr{'A'};
    auto v = view::repeat(chr);

    EXPECT_TRUE((std::ranges::Range<decltype(v)>));
    EXPECT_TRUE((std::ranges::InputRange<decltype(v)>));
    EXPECT_TRUE((std::ranges::ForwardRange<decltype(v)>));
    EXPECT_TRUE((std::ranges::BidirectionalRange<decltype(v)>));
    EXPECT_TRUE((std::ranges::RandomAccessRange<decltype(v)>));
    EXPECT_FALSE((std::ranges::ContiguousRange<decltype(v)>));
    EXPECT_TRUE((std::ranges::View<decltype(v)>));
    EXPECT_FALSE((std::ranges::SizedRange<decltype(v)>));
    EXPECT_FALSE((std::ranges::CommonRange<decltype(v)>));
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v), char>));
}

TEST(general, iterator)
{
    auto v = view::repeat('A');

    EXPECT_TRUE(v.begin() == v.begin());

    EXPECT_FALSE(v.begin() == v.end());
    EXPECT_FALSE(v.end() == v.begin());
    EXPECT_TRUE(v.begin() != v.end());
    EXPECT_TRUE(v.end() != v.begin());

    EXPECT_FALSE(v.cbegin() == v.cend());
    EXPECT_FALSE(v.cend() == v.cbegin());
    EXPECT_TRUE(v.cbegin() != v.cend());
    EXPECT_TRUE(v.cend() != v.cbegin());

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

    auto cit = v.cbegin();
    cit = v.begin(); // assign it from const it
}

TEST(general, subscript_operator)
{
    auto v = view::repeat('A');

    EXPECT_EQ(v[0], 'A');
    EXPECT_EQ(v[126], 'A');
    EXPECT_EQ(v[78634126], 'A');

    v[234] = 'X';

    EXPECT_EQ(v[0], 'X');
    EXPECT_EQ(v[126], 'X');
    EXPECT_EQ(v[78634126], 'X');
}

TEST(view, factory)
{
    // const char
    {
        char const chr{'X'};
        auto v = view::repeat(chr);
        EXPECT_EQ(*v.begin(), chr);
    }

    // string
    {
        std::string str{"foobar"};
        auto v = view::repeat(str);
        EXPECT_EQ(*v.begin(), str);
        EXPECT_EQ(v[2345], str);
    }

    // view
    {
        auto view = std::string{"foobar"} | view::persist | std::view::take(3);
        auto v = view::repeat(view);
        EXPECT_TRUE(std::ranges::equal(*v.begin(), view));
    }

    // combinability
    {
        std::string str{"foobar"};
        auto v = view::repeat(str) | view::take_exactly(3);
        EXPECT_EQ(*v.begin(), str);
        EXPECT_EQ(std::ranges::size(v), 3u);
    }
}

constexpr char constexpr_class_and_iterator()
{
    auto v = view::repeat('A');

    auto it = v.begin();
    ++it;
    it[234] = 'X';

    return *it;
}

constexpr char constexpr_view()
{
    char chr{'A'};
    auto v = view::repeat(chr);
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
