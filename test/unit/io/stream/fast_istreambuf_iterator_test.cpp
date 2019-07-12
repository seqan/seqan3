// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/std/iterator>

using namespace seqan3;

TEST(fast_istreambuf_iterator, concept)
{
    EXPECT_TRUE(std::InputIterator<detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::ForwardIterator<detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::BidirectionalIterator<detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::RandomAccessIterator<detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE((std::OutputIterator<detail::fast_istreambuf_iterator<char>, char>));
}

TEST(fast_istreambuf_iterator, construction)
{
    using type = detail::fast_istreambuf_iterator<char>;
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
    detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_EQ(*it, 't');
    ++it;
    EXPECT_EQ(*it, 'e');
    it++;
    EXPECT_EQ(*it, 's');
}

TEST(fast_istreambuf_iterator, comparison)
{
    std::istringstream str{"test\n"};
    detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_FALSE(it == std::ranges::default_sentinel);
    EXPECT_FALSE(std::ranges::default_sentinel == it);
    EXPECT_TRUE(it != std::ranges::default_sentinel);
    EXPECT_TRUE(std::ranges::default_sentinel != it);
}
