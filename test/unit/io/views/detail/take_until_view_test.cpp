// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <span>

#include <seqan3/io/views/detail/take_until_view.hpp>
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
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), char>)); // lost by single_pass_input

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
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    do_test(seqan3::detail::take_until, is_newline, "foo\nbar");
}

TEST(view_take_until, functor_fail)
{
    using namespace std::literals;

    std::string vec{"foo"};
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    EXPECT_RANGE_EQ("foo"sv, vec | seqan3::detail::take_until(is_newline));
}

TEST(view_take_until, concepts)
{
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    auto adapt = seqan3::detail::take_until(is_newline);
    do_concepts(adapt, true);

    // mutable adapters make the view lose const-iterability
    char c2 = 'a';
    auto adapt2 = seqan3::detail::take_until(
        [c2](char c) mutable
        {
            c2 = 'b';
            return c == '\n';
        });
    do_concepts(adapt2, false);
}

// ============================================================================
//  view_take_until_or_throw
// ============================================================================

TEST(view_take_until_or_throw, unix_eol)
{
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    do_test(seqan3::detail::take_until_or_throw, is_newline, "foo\nbar");
}

TEST(view_take_until_or_throw, functor_fail)
{
    std::string vec{"foo"};

    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    EXPECT_THROW(std::ranges::for_each(vec | seqan3::detail::take_until_or_throw(is_newline), [](auto &&) {}),
                 seqan3::unexpected_end_of_input);
}

TEST(view_take_until_or_throw, concepts)
{
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    do_concepts(seqan3::detail::take_until_or_throw(is_newline), true);
}

// ============================================================================
//  take_until_and_consume
// ============================================================================

TEST(take_until_and_consume, unix_eol)
{
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    do_test(seqan3::detail::take_until_and_consume, is_newline, "foo\n\n\n\nbar");
}

TEST(take_until_and_consume, consume)
{
    using namespace std::literals;

    std::string vec{"foo\n\n\n\nbar"};
    auto input_view = vec | seqan3::views::single_pass_input;

    auto take_until = input_view
                    | seqan3::detail::take_until_and_consume(
                          [](char c)
                          {
                              return c == '\n';
                          });

    // consumes "foo\n\n\n\n"
    EXPECT_RANGE_EQ("foo"sv, take_until);

    // next char in input range should be 'b'
    EXPECT_EQ(*input_view.begin(), 'b');
}

TEST(take_until_and_consume, functor_fail)
{
    using namespace std::literals;

    std::string vec{"foo"};
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    EXPECT_RANGE_EQ("foo"sv, vec | seqan3::detail::take_until_and_consume(is_newline));
}

TEST(take_until_and_consume, concepts)
{
    auto is_newline = [](char c)
    {
        return c == '\n';
    };
    do_concepts(seqan3::detail::take_until_and_consume(is_newline), true);
}
