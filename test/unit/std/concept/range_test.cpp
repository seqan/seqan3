// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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
    EXPECT_TRUE ((std::ranges::Range<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::Range<std::list<char>>));
    EXPECT_TRUE ((std::ranges::Range<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::Range<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::Range<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::Range<std::string>));

    EXPECT_TRUE ((std::ranges::Range<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::Range<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::Range<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::Range<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::Range<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::Range<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::Range<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::Range<std::vector<char> const &>));
}

TEST(range_concepts, InputRange)
{
    EXPECT_TRUE ((std::ranges::InputRange<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::InputRange<std::list<char>>));
    EXPECT_TRUE ((std::ranges::InputRange<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::InputRange<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::InputRange<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::InputRange<std::string>));

    EXPECT_TRUE ((std::ranges::InputRange<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::InputRange<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::InputRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::InputRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::InputRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::InputRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::InputRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::InputRange<std::vector<char> const &>));
}

TEST(range_concepts, ForwardRange)
{
    EXPECT_TRUE ((std::ranges::ForwardRange<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<std::list<char>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<std::string>));

    EXPECT_TRUE ((std::ranges::ForwardRange<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::ForwardRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::ForwardRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::ForwardRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::ForwardRange<std::vector<char> const &>));
}

TEST(range_concepts, BidirectionalRange)
{
    EXPECT_FALSE((std::ranges::BidirectionalRange<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::list<char>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::string>));

    EXPECT_TRUE ((std::ranges::BidirectionalRange<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::BidirectionalRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::BidirectionalRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::BidirectionalRange<std::vector<char> const &>));
}

TEST(range_concepts, RandomAccessRange)
{
    EXPECT_FALSE((std::ranges::RandomAccessRange<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::RandomAccessRange<std::list<char>>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<std::string>));

    EXPECT_TRUE ((std::ranges::RandomAccessRange<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::RandomAccessRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::RandomAccessRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((std::ranges::RandomAccessRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::RandomAccessRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::RandomAccessRange<std::vector<char> const &>));
}

TEST(range_concepts, ContiguousRange)
{
    EXPECT_FALSE((std::ranges::ContiguousRange<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::ContiguousRange<std::list<char>>));
    EXPECT_TRUE ((std::ranges::ContiguousRange<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::ContiguousRange<std::vector<char>>));
    EXPECT_FALSE((std::ranges::ContiguousRange<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::ContiguousRange<std::string>));

//     EXPECT_TRUE ((std::ranges::ContiguousRange<sdsl::int_vector<8>>));
    EXPECT_FALSE((std::ranges::ContiguousRange<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::ContiguousRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::ContiguousRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((std::ranges::ContiguousRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_FALSE((std::ranges::ContiguousRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::ContiguousRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::ContiguousRange<std::vector<char> const &>));
}

TEST(range_concepts, OutputRange)
{
    EXPECT_TRUE ((std::ranges::OutputRange<std::forward_list<char>, char>));
    EXPECT_TRUE ((std::ranges::OutputRange<std::list<char>, char>));
    EXPECT_TRUE ((std::ranges::OutputRange<std::array<char, 2>, char>));
    EXPECT_TRUE ((std::ranges::OutputRange<std::vector<char>, char>));
    EXPECT_TRUE ((std::ranges::OutputRange<std::deque<char>, char>));
    EXPECT_TRUE ((std::ranges::OutputRange<std::string, char>));

    EXPECT_TRUE ((std::ranges::OutputRange<sdsl::int_vector<8>, uint8_t>));
//     EXPECT_TRUE ((std::ranges::OutputRange<sdsl::int_vector<9>, uint8_t>));

    EXPECT_FALSE((std::ranges::OutputRange<ranges::any_view<char, ranges::category::input>, char>));
    EXPECT_FALSE((std::ranges::OutputRange<ranges::any_view<char, ranges::category::forward>, char>));
    EXPECT_FALSE((std::ranges::OutputRange<ranges::any_view<char, ranges::category::bidirectional>, char>));
    EXPECT_FALSE((std::ranges::OutputRange<ranges::any_view<char, ranges::category::random_access>, char>));

    EXPECT_TRUE ((std::ranges::OutputRange<std::vector<char> &, char>));
    EXPECT_FALSE((std::ranges::OutputRange<std::vector<char> const, char>));
    EXPECT_FALSE((std::ranges::OutputRange<std::vector<char> const &, char>));
}

TEST(range_concepts, SizedRange)
{
    EXPECT_FALSE((std::ranges::SizedRange<std::forward_list<char>>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::list<char>>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::array<char, 2>>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::vector<char>>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::deque<char>>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::string>));

    EXPECT_TRUE ((std::ranges::SizedRange<sdsl::int_vector<8>>));
    EXPECT_TRUE ((std::ranges::SizedRange<sdsl::int_vector<9>>));

    EXPECT_FALSE((std::ranges::SizedRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_FALSE((std::ranges::SizedRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_FALSE((std::ranges::SizedRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_FALSE((std::ranges::SizedRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::SizedRange<std::vector<char> &>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::SizedRange<std::vector<char> const &>));
}

TEST(range_concepts, View)
{
    EXPECT_FALSE((std::ranges::View<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::View<std::list<char>>));
    EXPECT_FALSE((std::ranges::View<std::array<char, 2>>));
    EXPECT_FALSE((std::ranges::View<std::vector<char>>));
    EXPECT_FALSE((std::ranges::View<std::deque<char>>));
    EXPECT_FALSE((std::ranges::View<std::string>));

    EXPECT_FALSE((std::ranges::View<sdsl::int_vector<8>>));
    EXPECT_FALSE((std::ranges::View<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::View<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::View<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::View<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::View<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_FALSE((std::ranges::View<std::vector<char> &>));
    EXPECT_FALSE((std::ranges::View<std::vector<char> const>));
    EXPECT_FALSE((std::ranges::View<std::vector<char> const &>));
}

TEST(range_concepts, ViewableRange)
{
    EXPECT_FALSE((std::ranges::ViewableRange<std::forward_list<char>>));
    EXPECT_FALSE((std::ranges::ViewableRange<std::list<char>>));
    EXPECT_FALSE((std::ranges::ViewableRange<std::array<char, 2>>));
    EXPECT_FALSE((std::ranges::ViewableRange<std::vector<char>>));
    EXPECT_FALSE((std::ranges::ViewableRange<std::deque<char>>));
    EXPECT_FALSE((std::ranges::ViewableRange<std::string>));

    EXPECT_FALSE((std::ranges::ViewableRange<sdsl::int_vector<8>>));
    EXPECT_FALSE((std::ranges::ViewableRange<sdsl::int_vector<9>>));

    EXPECT_TRUE ((std::ranges::ViewableRange<ranges::any_view<char, ranges::category::input>>));
    EXPECT_TRUE ((std::ranges::ViewableRange<ranges::any_view<char, ranges::category::forward>>));
    EXPECT_TRUE ((std::ranges::ViewableRange<ranges::any_view<char, ranges::category::bidirectional>>));
    EXPECT_TRUE ((std::ranges::ViewableRange<ranges::any_view<char, ranges::category::random_access>>));

    EXPECT_TRUE ((std::ranges::ViewableRange<std::vector<char> &>));
    EXPECT_FALSE((std::ranges::ViewableRange<std::vector<char> const>));
    EXPECT_TRUE ((std::ranges::ViewableRange<std::vector<char> const &>));
}
