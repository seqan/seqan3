// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
