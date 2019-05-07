// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/view/unique.hpp>

#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t, typename fun_t>
void do_test(adaptor_t const & adaptor, fun_t && fun, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor(fun);
    EXPECT_EQ("foo", std::string(v));

    // function notation
    std::string v2 = adaptor(vec, fun);
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | adaptor(fun) | ranges::view::unique;
    EXPECT_EQ("fo", std::string(v3));
    std::string v3b = vec | std::view::reverse | adaptor(fun) | ranges::view::unique;
    EXPECT_EQ("rab", v3b);

    // pointer as iterator
    std::span s{std::ranges::data(vec), vec.size()};
    auto v4 = s | adaptor(fun);
    EXPECT_EQ("foo", std::string(v4));

    // comparability against self
    EXPECT_TRUE(std::ranges::equal(v,v));
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const_it)
{
    std::string vec{"foo\nbar"};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(ConstIterableRange<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), char>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_FALSE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_EQ(ConstIterableRange<decltype(v1)>, const_it);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), char>));

    auto v2 = vec | view::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::ForwardRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::BidirectionalRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::RandomAccessRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::View<decltype(v2)>);
    EXPECT_FALSE(std::ranges::SizedRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v2)>);
    EXPECT_FALSE(ConstIterableRange<decltype(v2)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v2), char>));
}

// ============================================================================
//  view_take_until
// ============================================================================

TEST(view_take_until, unix_eol)
{
    do_test(view::take_until, [] (char c) { return c == '\n'; }, "foo\nbar");
}

TEST(view_take_until, functor_fail)
{
    std::string vec{"foo"};
    std::string v;
    EXPECT_NO_THROW(( v = vec | view::take_until([] (char c) { return c == '\n'; }) ));
    EXPECT_EQ("foo", v);
}

TEST(view_take_until, concepts)
{
    auto adapt = view::take_until([] (char c) { return c == '\n'; });
    do_concepts(adapt, true);

    // mutable adapters make the view loose const-iterability, but this is not checked by conepts unfortunately
//     auto adapt2 = view::take_until([count = 0] (char c) mutable { ++count; return c == '\n'; });
//     do_concepts(adapt2, false);
}

// ============================================================================
//  view_take_until_or_throw
// ============================================================================

TEST(view_take_until_or_throw, unix_eol)
{
    do_test(view::take_until_or_throw, [] (char c) { return c == '\n'; }, "foo\nbar");
}

TEST(view_take_until_or_throw, functor_fail)
{
    std::string vec{"foo"};
    EXPECT_THROW(std::string v = vec | view::take_until_or_throw([] (char c) { return c == '\n'; }),
                 unexpected_end_of_input);
}

TEST(view_take_until_or_throw, concepts)
{
    do_concepts(view::take_until_or_throw([] (char c) { return c == '\n'; }), true);
}
