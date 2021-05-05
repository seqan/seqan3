// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

#include <seqan3/io/detail/take_until_view.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t, typename fun_t>
void do_test(adaptor_t const & adaptor, fun_t && fun, std::string const & vec)
{
    using namespace std::literals;

    // pipe notation
    EXPECT_RANGE_EQ("foo"sv, vec | adaptor(fun));

    // function notation
    EXPECT_RANGE_EQ("foo"sv, adaptor(vec, fun));

    // combinability
    EXPECT_RANGE_EQ("fo"sv, vec | adaptor(fun) | std::views::take(2));
    EXPECT_RANGE_EQ("rab"sv, vec | std::views::reverse | adaptor(fun) | std::views::take(3));

    // Test combinability of take_until with std::reverse as second view, this caused a problem here:
    // https://github.com/seqan/seqan3/issues/1754
    EXPECT_RANGE_EQ("oof"sv, vec | adaptor(fun) | std::views::reverse);

    // pointer as iterator
    std::span s{std::ranges::data(vec), vec.size()};
    EXPECT_RANGE_EQ("foo"sv, s | adaptor(fun));
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
    auto is_newline = [] (char c) { return c == '\n'; };
    do_test(seqan3::detail::take_until, is_newline, "foo\nbar");
}

TEST(view_take_until, functor_fail)
{
    using namespace std::literals;

    std::string vec{"foo"};
    auto is_newline = [](char c){return c == '\n'; };
    EXPECT_RANGE_EQ("foo"sv, vec | seqan3::detail::take_until(is_newline));
}

TEST(view_take_until, concepts)
{
    auto is_newline = [] (char c) { return c == '\n'; };
    auto adapt = seqan3::detail::take_until(is_newline);
    do_concepts(adapt, true);

    // mutable adapters make the view loose const-iterability, but this is not checked by conepts unfortunately
//     auto adapt2 = seqan3::detail::take_until([count = 0] (char c) mutable { ++count; return c == '\n'; });
//     do_concepts(adapt2, false);
}

// ============================================================================
//  view_take_until_or_throw
// ============================================================================

TEST(view_take_until_or_throw, unix_eol)
{
    auto is_newline = [] (char c) { return c == '\n'; };
    do_test(seqan3::detail::take_until_or_throw, is_newline, "foo\nbar");
}

TEST(view_take_until_or_throw, functor_fail)
{
    std::string vec{"foo"};

    auto is_newline = [] (char c) { return c == '\n'; };
    EXPECT_THROW(std::ranges::for_each(vec | seqan3::detail::take_until_or_throw(is_newline), [](auto &&){}),
                 seqan3::unexpected_end_of_input);
}

TEST(view_take_until_or_throw, concepts)
{
    auto is_newline = [](char c){return c == '\n'; };
    do_concepts(seqan3::detail::take_until_or_throw(is_newline), true);
}
