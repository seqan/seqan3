// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <array>
#include <deque>
#include <forward_list>
#include <list>
#include <string>
#include <vector>

#include <range/v3/view/any_view.hpp>

#include <sdsl/int_vector.hpp>

#include <seqan3/std/ranges>

TEST(range_concepts, Range)
{
    EXPECT_TRUE ((std::ranges::range<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::range<std::string>));

    EXPECT_TRUE ((std::ranges::range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::range<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::range<std::vector<char> const &>));
}

TEST(range_concepts, input_range)
{
    EXPECT_TRUE ((std::ranges::input_range<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::input_range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::input_range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::input_range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::input_range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::input_range<std::string>));

    EXPECT_TRUE ((std::ranges::input_range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::input_range<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::input_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::input_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::input_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::input_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::input_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::input_range<std::vector<char> const &>));
}

TEST(range_concepts, forward_range)
{
    EXPECT_TRUE ((std::ranges::forward_range<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::forward_range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::forward_range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::forward_range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::forward_range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::forward_range<std::string>));

    EXPECT_TRUE ((std::ranges::forward_range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::forward_range<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::forward_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::forward_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::forward_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::forward_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::forward_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::forward_range<std::vector<char> const &>));
}

TEST(range_concepts, bidirectional_range)
{
    EXPECT_FALSE((std::ranges::bidirectional_range<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<std::string>));

    EXPECT_TRUE ((std::ranges::bidirectional_range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::bidirectional_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::bidirectional_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::bidirectional_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::bidirectional_range<std::vector<char> const &>));
}

TEST(range_concepts, random_access_range)
{
    EXPECT_FALSE((std::ranges::random_access_range<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::random_access_range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::random_access_range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::random_access_range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::random_access_range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::random_access_range<std::string>));

    EXPECT_TRUE ((std::ranges::random_access_range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::random_access_range<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::random_access_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::random_access_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((std::ranges::random_access_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::random_access_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::random_access_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::random_access_range<std::vector<char> const &>));
}

TEST(range_concepts, contiguous_range)
{
    EXPECT_FALSE((std::ranges::contiguous_range<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::contiguous_range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::contiguous_range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::contiguous_range<std::vector<char>>));
    EXPECT_FALSE((std::ranges::contiguous_range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::contiguous_range<std::string>));

//     EXPECT_TRUE ((std::ranges::contiguous_range<sdsl::int_vector<8>>));
    EXPECT_FALSE((std::ranges::contiguous_range<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::contiguous_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::contiguous_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((std::ranges::contiguous_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_FALSE((std::ranges::contiguous_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::contiguous_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::contiguous_range<std::vector<char> const &>));
}

TEST(range_concepts, output_range)
{
    EXPECT_TRUE ((std::ranges::output_range<std::forward_list<char>, char>));
    EXPECT_TRUE ((std::ranges::output_range<std::list<char>, char>));
    EXPECT_TRUE ((std::ranges::output_range<std::array<char, 2>, char>));
    EXPECT_TRUE ((std::ranges::output_range<std::vector<char>, char>));
    EXPECT_TRUE ((std::ranges::output_range<std::deque<char>, char>));
    EXPECT_TRUE ((std::ranges::output_range<std::string, char>));

    EXPECT_TRUE ((std::ranges::output_range<sdsl::int_vector<8>, uint8_t>));
//     EXPECT_TRUE ((std::ranges::output_range<sdsl::int_vector<9>, uint8_t>));

    EXPECT_FALSE((std::ranges::output_range<ranges::any_view<char, ranges::category::input>, char>));
    EXPECT_FALSE((std::ranges::output_range<ranges::any_view<char, ranges::category::forward>, char>));
    EXPECT_FALSE((std::ranges::output_range<ranges::any_view<char, ranges::category::bidirectional>, char>));
    EXPECT_FALSE((std::ranges::output_range<ranges::any_view<char, ranges::category::random_access>, char>));

    EXPECT_TRUE ((std::ranges::output_range<std::vector<char> &, char>));
    EXPECT_FALSE((std::ranges::output_range<std::vector<char> const, char>));
    EXPECT_FALSE((std::ranges::output_range<std::vector<char> const &, char>));
}

TEST(range_concepts, sized_range)
{
    EXPECT_FALSE((std::ranges::sized_range<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::sized_range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::sized_range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::sized_range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::sized_range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::sized_range<std::string>));

    EXPECT_TRUE ((std::ranges::sized_range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::sized_range<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::sized_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::sized_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((std::ranges::sized_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_FALSE((std::ranges::sized_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::sized_range<std::vector<char> &>));
    EXPECT_TRUE ((std::ranges::sized_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::sized_range<std::vector<char> const &>));
}

TEST(range_concepts, View)
{
    EXPECT_FALSE((std::ranges::view<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::view<std::list<char>>));
    EXPECT_FALSE((std::ranges::view<std::array<char, 2>>));
    EXPECT_FALSE((std::ranges::view<std::vector<char>>));
    EXPECT_FALSE((std::ranges::view<std::deque<char>>));
    EXPECT_FALSE((std::ranges::view<std::string>));

    EXPECT_FALSE((std::ranges::view<sdsl::int_vector<8>>));
    EXPECT_FALSE((std::ranges::view<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::view<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::view<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::view<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::view<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_FALSE((std::ranges::view<std::vector<char> &>));
    EXPECT_FALSE((std::ranges::view<std::vector<char> const>));
    EXPECT_FALSE((std::ranges::view<std::vector<char> const &>));
}

TEST(range_concepts, viewable_range)
{
    EXPECT_FALSE((std::ranges::viewable_range<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::viewable_range<std::list<char>>));
    EXPECT_FALSE((std::ranges::viewable_range<std::array<char, 2>>));
    EXPECT_FALSE((std::ranges::viewable_range<std::vector<char>>));
    EXPECT_FALSE((std::ranges::viewable_range<std::deque<char>>));
    EXPECT_FALSE((std::ranges::viewable_range<std::string>));

    EXPECT_FALSE((std::ranges::viewable_range<sdsl::int_vector<8>>));
    EXPECT_FALSE((std::ranges::viewable_range<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::viewable_range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::viewable_range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::viewable_range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::viewable_range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::viewable_range<std::vector<char> &>));
    EXPECT_FALSE((std::ranges::viewable_range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::viewable_range<std::vector<char> const &>));
}
