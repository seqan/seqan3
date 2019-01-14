// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/view/reverse.hpp>
#include <seqan3/std/view/common.hpp>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor(3);
    EXPECT_EQ("foo", std::string(v));

    // function notation
    std::string v2 = adaptor(vec, 3);
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | adaptor(3) | adaptor(3) | ranges::view::unique;
    EXPECT_EQ("fo", std::string(v3));
    std::string v3b = vec | view::reverse | adaptor(3) | ranges::view::unique;
    EXPECT_EQ("rab", v3b);
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const exactly)
{
    std::string vec{"foobar"};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), char>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(v1)>, exactly);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), char>));

    auto v2 = vec | view::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::ForwardRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::BidirectionalRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::RandomAccessRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::View<decltype(v2)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(v2)>, exactly);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v2)>);
    EXPECT_FALSE(const_iterable_concept<decltype(v2)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v2), char>));
}

// ============================================================================
//  view_take
// ============================================================================

TEST(view_take, regular)
{
    do_test(view::take, "foobar");
}

TEST(view_take, concepts)
{
    do_concepts(view::take(3), false);
}

TEST(view_take, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_NO_THROW(( view::take(vec, 4) )); // no parsing

    std::string v;
    EXPECT_NO_THROW(( v = vec | view::single_pass_input | view::take(4) )); // full parsing on conversion
    EXPECT_EQ("foo", v);
}

// ============================================================================
//  view_take_exactly
// ============================================================================

TEST(view_take_exactly, regular)
{
    do_test(view::take_exactly, "foobar");
}

TEST(view_take_exactly, concepts)
{
    do_concepts(view::take_exactly(3), true);
}

TEST(view_take_exactly, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_NO_THROW(( view::take_exactly(vec, 4) )); // no parsing

    std::string v;
    EXPECT_NO_THROW(( v = vec | view::single_pass_input | view::take_exactly(4) )); // full parsing on conversion
    EXPECT_EQ("foo", v);

    auto v2 = vec | view::single_pass_input | view::take_exactly(4);
    EXPECT_EQ(size(v2), 4u); // here be dragons
}

// ============================================================================
//  view_take_exactly_or_throw
// ============================================================================

TEST(view_take_exactly_or_throw, regular)
{
    do_test(view::take_exactly_or_throw, "foo\nbar");
}

TEST(view_take_exactly_or_throw, concepts)
{
    do_concepts(view::take_exactly_or_throw(3), true);
}

TEST(view_take_exactly_or_throw, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_THROW(( view::take_exactly_or_throw(vec, 4) ),
                   std::invalid_argument); // no parsing, but throws on construction

    std::string v;
    EXPECT_THROW(( v = vec | view::single_pass_input | view::take_exactly_or_throw(4)),
                   unexpected_end_of_input); // full parsing on conversion, throw on conversion
}
