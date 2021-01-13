// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/iterator>

#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>

TEST(fast_istreambuf_iterator, concept)
{
    EXPECT_TRUE(std::input_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::forward_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::bidirectional_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::random_access_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE((std::output_iterator<seqan3::detail::fast_istreambuf_iterator<char>, char>));
}

TEST(fast_istreambuf_iterator, construction)
{
    using type = seqan3::detail::fast_istreambuf_iterator<char>;
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<type>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<type>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<type>);
    EXPECT_TRUE((std::is_constructible_v<type, std::basic_streambuf<char> &>));
}

TEST(fast_istreambuf_iterator, basic)
{
    std::istringstream str{"test"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_EQ(*it, 't');
    ++it;
    EXPECT_EQ(*it, 'e');
    it++;
    EXPECT_EQ(*it, 's');
}

TEST(fast_istreambuf_iterator, comparison)
{
    std::istringstream str{"test\n"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_FALSE(it == std::default_sentinel);
    EXPECT_FALSE(std::default_sentinel == it);
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_TRUE(std::default_sentinel != it);
}
