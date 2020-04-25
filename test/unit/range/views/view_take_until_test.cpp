// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/view/unique.hpp>

#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/std/span>
#include <seqan3/test/expect_range_eq.hpp>

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t, typename fun_t>
void do_test(adaptor_t const & adaptor, fun_t && fun, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor(fun);
    EXPECT_EQ("foo", v  | seqan3::views::to<std::string>);

    // function notation
    std::string v2 = adaptor(vec, fun) | seqan3::views::to<std::string>;
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | adaptor(fun) | ranges::view::unique;
    EXPECT_EQ("fo", v3  | seqan3::views::to<std::string>);
    std::string v3b = vec | std::views::reverse | adaptor(fun) | ranges::view::unique | seqan3::views::to<std::string>;
    EXPECT_EQ("rab", v3b);

    // pointer as iterator
    std::span s{std::ranges::data(vec), vec.size()};
    auto v4 = s | adaptor(fun);
    EXPECT_EQ("foo", v4 | seqan3::views::to<std::string>);

    // comparability against self
    EXPECT_RANGE_EQ(v, v);
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const_it)
{
    std::string vec{"foo\nbar"};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), char>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_EQ(seqan3::const_iterable_range<decltype(v1)>, const_it);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));

    auto v2 = vec | seqan3::views::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v2), char>));

    // explicit test for non const-iterable views
    // https://github.com/seqan/seqan3/pull/1734#discussion_r408829267
    auto const & v2_cref = v2;

    EXPECT_FALSE(std::ranges::input_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::view<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2_cref)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2_cref)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2_cref), int>));
}

// ============================================================================
//  view_take_until
// ============================================================================

TEST(view_take_until, unix_eol)
{
    do_test(seqan3::views::take_until, [] (char c) { return c == '\n'; }, "foo\nbar");
}

TEST(view_take_until, functor_fail)
{
    std::string vec{"foo"};
    std::string v;
    EXPECT_NO_THROW(( v = vec
                        | seqan3::views::take_until([] (char c) { return c == '\n'; })
                        | seqan3::views::to<std::string> ));
    EXPECT_EQ("foo", v);
}

TEST(view_take_until, concepts)
{
    auto adapt = seqan3::views::take_until([] (char c) { return c == '\n'; });
    do_concepts(adapt, true);

    // mutable adapters make the view loose const-iterability, but this is not checked by conepts unfortunately
//     auto adapt2 = seqan3::views::take_until([count = 0] (char c) mutable { ++count; return c == '\n'; });
//     do_concepts(adapt2, false);
}

// ============================================================================
//  view_take_until_or_throw
// ============================================================================

TEST(view_take_until_or_throw, unix_eol)
{
    do_test(seqan3::views::take_until_or_throw, [] (char c) { return c == '\n'; }, "foo\nbar");
}

TEST(view_take_until_or_throw, functor_fail)
{
    std::string vec{"foo"};
    EXPECT_THROW(std::string v = vec
                               | seqan3::views::take_until_or_throw([] (char c) { return c == '\n'; })
                               | seqan3::views::to<std::string>,
                 seqan3::unexpected_end_of_input);
}

TEST(view_take_until_or_throw, concepts)
{
    do_concepts(seqan3::views::take_until_or_throw([] (char c) { return c == '\n'; }), true);
}
