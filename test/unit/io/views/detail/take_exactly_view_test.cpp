// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <deque>
#include <list>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

#include "../../../range/iterator_test_template.hpp"
#include "../../../range/range_test_template.hpp"

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    using namespace std::literals;

    // pipe notation
    EXPECT_RANGE_EQ("foo"sv, vec | adaptor(3));

    // function notation
    EXPECT_RANGE_EQ("foo"sv, adaptor(vec, 3));

    // combinability
    EXPECT_RANGE_EQ("fo"sv, vec | adaptor(3) | adaptor(3) | std::views::take(2));
    EXPECT_RANGE_EQ("rab"sv, vec | std::views::reverse | adaptor(3) | std::views::take(3));
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const exactly)
{
    std::vector vec{1, 2, 3};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::contiguous_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), int>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::contiguous_range<decltype(v1)>); // same iterator/sentinel as from vector
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), int>));

    auto v3 = vec
            | std::views::transform(
                  [](auto && v)
                  {
                      return v;
                  })
            | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v3)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::view<decltype(v3)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v3)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v3)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v3), int>));

    auto v2 = vec | seqan3::views::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(v2)>, exactly);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), int>));

    // explicit test for non const-iterable views
    // https://github.com/seqan/seqan3/pull/1734#discussion_r408829267
    auto const & v2_cref = v2;

    EXPECT_FALSE(std::ranges::input_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::view<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2_cref)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2_cref)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2_cref), int>));
}

// ============================================================================
//  view_take_exactly
// ============================================================================

TEST(view_take_exactly, regular)
{
    do_test(seqan3::detail::take_exactly, "foobar");
}

TEST(view_take_exactly, concepts)
{
    do_concepts(seqan3::detail::take_exactly(3), true);
}

TEST(view_take_exactly, underlying_is_shorter)
{
    using namespace std::literals;

    std::string vec{"foo"};
    EXPECT_NO_THROW((seqan3::detail::take_exactly(vec, 4))); // no parsing

    // full parsing on conversion
    EXPECT_RANGE_EQ("foo"sv, vec | seqan3::views::single_pass_input | seqan3::detail::take_exactly(4));

    auto v2 = vec | seqan3::views::single_pass_input | seqan3::detail::take_exactly(4);
    EXPECT_EQ(std::ranges::size(v2), 4u); // here be dragons
}

TEST(view_take_exactly, shrink_size_on_input_ranges)
{
    std::string vec{"foobar"};
    auto v = vec | seqan3::views::single_pass_input | seqan3::detail::take_exactly(3);

    EXPECT_EQ(std::ranges::size(v), 3u);
    EXPECT_EQ(*std::ranges::begin(v), 'f');

    auto it = std::ranges::begin(v);
    ++it;

    EXPECT_EQ(std::ranges::size(v), 2u);
    EXPECT_EQ(*std::ranges::begin(v), 'o');

    ++it;
    ++it;

    EXPECT_EQ(std::ranges::size(v), 0u); // view is empty now
}

struct view_take_exactly1_test_fixture : public range_test_fixture
{
    using range_value_t = char;
    using range_reference_t = char const &;

    using range_const_value_t = char;
    using range_const_reference_t = char const &;

    static constexpr bool input_range = true;
    static constexpr bool forward_range = true;
    static constexpr bool bidirectional_range = true;
    static constexpr bool random_access_range = true;
    static constexpr bool contiguous_range = true;
    static constexpr bool output_range = false;

    static constexpr bool common_range = true;
    static constexpr bool viewable_range = true;
    static constexpr bool view = true;
    static constexpr bool sized_range = true;
    static constexpr bool const_iterable_range = true;
    static constexpr bool size_member = true;
    static constexpr bool const_size_member = true;
    static constexpr bool subscript_member = true;

    std::string_view expected_range()
    {
        return {"01234"};
    }

    auto range()
    {
        return std::string_view{"0123456789"} | seqan3::detail::take_exactly(5);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(view_take_exactly1_test, range_test, view_take_exactly1_test_fixture, );
INSTANTIATE_TYPED_TEST_SUITE_P(view_take_exactly1_test, iterator_fixture, view_take_exactly1_test_fixture, );

// This tests a use cases from format_bam_test where std::ranges::subrange assumed const_iterable even though
// seqan3::detail::take_exactly SHOULD lose this property on input_iterator (in this case it manages the data within the
// view and needs mutable access).
struct view_take_exactly2_test_fixture : public range_test_fixture
{
    using range_value_t = char;
    using range_reference_t = char;

    static constexpr bool input_range = true;
    static constexpr bool forward_range = false;
    static constexpr bool bidirectional_range = false;
    static constexpr bool random_access_range = false;
    static constexpr bool contiguous_range = false;
    static constexpr bool output_range = false;

    static constexpr bool common_range = false;
    static constexpr bool viewable_range = true;
    static constexpr bool view = true;
    static constexpr bool sized_range = true;           // seqan3::detail::take_exactly adds this property
    static constexpr bool const_iterable_range = false; // seqan3::detail::take_exactly loses this property
    static constexpr bool size_member = true;           // seqan3::detail::take_exactly adds this property
    static constexpr bool const_size_member = true;     // seqan3::detail::take_exactly adds this property
    static constexpr bool subscript_member = false;

    std::string_view expected_range()
    {
        return "01234";
    }

    std::string _range{"0123456789"};

    auto range()
    {
        static std::istringstream istream{};
        istream.str(_range); // reset underlying stream on each invocation

        using iterator_t = seqan3::detail::fast_istreambuf_iterator<char>;
        using sentinel_t = std::default_sentinel_t;
        using subrange_t = std::ranges::subrange<iterator_t, sentinel_t>;

        return subrange_t{iterator_t{*istream.rdbuf()}, sentinel_t{}} | seqan3::detail::take_exactly(5);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(view_take_exactly2_test, range_test, view_take_exactly2_test_fixture, );
INSTANTIATE_TYPED_TEST_SUITE_P(view_take_exactly2_test, iterator_fixture, view_take_exactly2_test_fixture, );

// ============================================================================
//  view_take_exactly_or_throw
// ============================================================================

TEST(view_take_exactly_or_throw, regular)
{
    do_test(seqan3::detail::take_exactly_or_throw, "foo\nbar");
}

TEST(view_take_exactly_or_throw, concepts)
{
    do_concepts(seqan3::detail::take_exactly_or_throw(3), true);
}

TEST(view_take_exactly_or_throw, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_THROW((seqan3::detail::take_exactly_or_throw(vec, 4)),
                 std::invalid_argument); // no parsing, but throws in adaptor

    std::list l{'f', 'o', 'o'};
    EXPECT_THROW((seqan3::detail::view_take_exactly<std::views::all_t<std::list<char> &>, true>(l, 4)),
                 std::invalid_argument); // no parsing, but throws on construction

    EXPECT_THROW(
        std::ranges::for_each(vec | seqan3::views::single_pass_input | seqan3::detail::take_exactly_or_throw(4),
                              [](auto &&) {}),
        seqan3::unexpected_end_of_input); // full parsing on conversion, throw on conversion
}
