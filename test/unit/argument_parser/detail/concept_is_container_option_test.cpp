// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/argument_parser/detail/concept.hpp>

TEST(is_container_option_concept_test, std_vector)
{
    EXPECT_TRUE(seqan3::detail::is_container_option<std::vector<int>>);
    EXPECT_TRUE(seqan3::detail::is_container_option<std::vector<char>>);
    EXPECT_TRUE(seqan3::detail::is_container_option<std::vector<int> &>);
    EXPECT_TRUE(seqan3::detail::is_container_option<std::vector<char> &>);
}

TEST(is_container_option_concept_test, std_string)
{
    EXPECT_FALSE(seqan3::detail::is_container_option<std::string>);
    EXPECT_FALSE(seqan3::detail::is_container_option<std::string &>);
}
