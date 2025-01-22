// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ranges>

#include <seqan3/io/views/detail/take_line_view.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    using namespace std::literals;

    // pipe notation
    EXPECT_RANGE_EQ("foo"sv, vec | adaptor);

    // function notation
    EXPECT_RANGE_EQ("foo"sv, adaptor(vec));

    // combinability
    EXPECT_RANGE_EQ("fo"sv, vec | adaptor | std::views::take(2));
    EXPECT_RANGE_EQ("rab"sv, vec | std::views::reverse | adaptor | std::views::take(3));

    // consuming behaviour
    auto v4 = vec | seqan3::views::single_pass_input;
    auto v5 = std::move(v4) | adaptor;
    EXPECT_RANGE_EQ("foo"sv, v5);
    EXPECT_EQ('b', *begin(v5)); // not newline
}

template <typename adaptor_t>
void do_concepts(adaptor_t const & adaptor)
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
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
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
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), char>));
}

// ============================================================================
//  view_take_line
// ============================================================================

TEST(view_take_line, unix_eol)
{
    do_test(seqan3::detail::take_line, "foo\nbar");
}

TEST(view_take_line, windows_eol)
{
    do_test(seqan3::detail::take_line, "foo\r\nbar");
}

TEST(view_take_line, no_eol)
{
    using namespace std::literals;

    std::string vec{"foo"};
    EXPECT_RANGE_EQ("foo"sv, vec | seqan3::detail::take_line);
}

TEST(view_take_line, eol_at_first_position)
{
    using namespace std::literals;

    using sbt = std::istreambuf_iterator<char>;

    std::istringstream vec{"\n\nfoo"};
    auto stream_view = std::ranges::subrange<decltype(sbt{vec}), decltype(sbt{})>{sbt{vec}, sbt{}};

    EXPECT_RANGE_EQ(""sv, stream_view | seqan3::detail::take_line);
    EXPECT_RANGE_EQ("foo"sv, stream_view | seqan3::detail::take_line);
}

TEST(view_take_line, concepts)
{
    do_concepts(seqan3::detail::take_line);
}

// ============================================================================
//  view_take_line_or_throw
// ============================================================================

TEST(view_take_line_or_throw, unix_eol)
{
    do_test(seqan3::detail::take_line_or_throw, "foo\nbar");
}

TEST(view_take_line_or_throw, windows_eol)
{
    do_test(seqan3::detail::take_line_or_throw, "foo\r\nbar");
}

TEST(view_take_line_or_throw, no_eol)
{
    std::string vec{"foo"};
    EXPECT_THROW(std::ranges::for_each(vec | seqan3::detail::take_line_or_throw, [](auto &&) {}),
                 seqan3::unexpected_end_of_input);
}

TEST(view_take_line_or_throw, concepts)
{
    do_concepts(seqan3::detail::take_line_or_throw);
}

// ============================================================================
//  bug
// ============================================================================

TEST(view_take_line, reverse_bug)
{
    using namespace std::literals;

    std::string vec{"foo\nbar"};
    auto v1 = vec | seqan3::detail::take_line;
    EXPECT_RANGE_EQ("foo"sv, v1);

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));

    // No build failure, but wrong results:
    //     auto v2 = v1 | std::views::reverse;
    //     EXPECT_EQ("oof", std::string(v2));
}
