// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/iterator>

#include <seqan3/io/detail/ignore_output_iterator.hpp>

using namespace std::literals;

TEST(ignore_output_iterator, concept_check)
{
    EXPECT_TRUE((std::output_iterator<seqan3::detail::ignore_output_iterator, char>));
    EXPECT_TRUE((std::output_iterator<seqan3::detail::ignore_output_iterator, int>));
    EXPECT_FALSE((std::input_iterator<seqan3::detail::ignore_output_iterator>));
}

TEST(ignore_output_iterator, assign)
{
    seqan3::detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(it = 'A'), std::addressof(it));
    EXPECT_EQ(std::addressof(it = 10), std::addressof(it));
}

TEST(ignore_output_iterator, pre_increment)
{
    seqan3::detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(++it), std::addressof(it));
}

TEST(ignore_output_iterator, post_increment)
{
    seqan3::detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(it++), std::addressof(it));
}

TEST(ignore_output_iterator, dereference)
{
    seqan3::detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(*it), std::addressof(it));
}
