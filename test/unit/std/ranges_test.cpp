// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/ranges>
#include <span>
#include <string>
#include <string_view>

#include <seqan3/test/expect_same_type.hpp>

TEST(ranges_test, take_view)
{
    std::string s{};

#if !SEQAN3_WORKAROUND_GCC_100139
    EXPECT_SAME_TYPE(decltype(std::views::take(std::span<int>{}, 0)), std::span<int>);
    EXPECT_SAME_TYPE(decltype(std::views::take(std::string_view{}, 0)), std::string_view);
    EXPECT_SAME_TYPE(decltype(std::views::take(std::views::empty<int>, 0)), std::ranges::empty_view<int>);
    EXPECT_SAME_TYPE(decltype(std::views::take(std::views::iota(0, 1), 0)), decltype(std::views::iota(0, 1)));

    EXPECT_SAME_TYPE(decltype(std::views::take(s, 0)),
                     (std::ranges::subrange<std::string::iterator, std::string::iterator>));
#endif // !SEQAN3_WORKAROUND_GCC_100139

    EXPECT_TRUE(std::ranges::borrowed_range<decltype(std::views::take(s, 0))>);
    EXPECT_TRUE(std::ranges::viewable_range<decltype(std::views::take(s, 0))>);
    EXPECT_TRUE(std::ranges::view<decltype(std::views::take(s, 0))>);
}

TEST(ranges_test, drop_view)
{
    std::string s{};

#if !SEQAN3_WORKAROUND_GCC_100139
    EXPECT_SAME_TYPE(decltype(std::views::drop(std::span<int>{}, 0)), std::span<int>);
    EXPECT_SAME_TYPE(decltype(std::views::drop(std::string_view{}, 0)), std::string_view);
    EXPECT_SAME_TYPE(decltype(std::views::drop(std::views::empty<int>, 0)), std::ranges::empty_view<int>);
    EXPECT_SAME_TYPE(decltype(std::views::drop(std::views::iota(0, 1), 0)), decltype(std::views::iota(0, 1)));

    EXPECT_SAME_TYPE(decltype(std::views::drop(s, 0)),
                     (std::ranges::subrange<std::string::iterator, std::string::iterator>));
#endif // !SEQAN3_WORKAROUND_GCC_100139

    EXPECT_TRUE(std::ranges::borrowed_range<decltype(std::views::drop(s, 0))>);
    EXPECT_TRUE(std::ranges::viewable_range<decltype(std::views::drop(s, 0))>);
    EXPECT_TRUE(std::ranges::view<decltype(std::views::drop(s, 0))>);
}

// https://github.com/ericniebler/range-v3/issues/1514
TEST(ranges_test, gcc10bug_rangev3_1514)
{
    {
        auto iota = std::views::iota(0, 5);
        EXPECT_EQ(*ranges::begin(iota), 0);
        EXPECT_EQ(*std::ranges::begin(iota), 0);
    }
    {
        // https://github.com/ericniebler/range-v3/issues/1514
        auto iota = std::views::iota(size_t{0u}, size_t{5u});
        EXPECT_EQ(*ranges::begin(iota), 0u);
        EXPECT_EQ(*std::ranges::begin(iota), 0u);
    }
}

// https://github.com/seqan/product_backlog/issues/372
TEST(ranges_test, issue372)
{
#ifdef __cpp_lib_ranges // >= gcc-10, range-v3 bug
    std::string vec{};
    std::istringstream istringstream{vec};

    auto v1 = std::ranges::subrange{std::istream_iterator<char>{istringstream}, std::default_sentinel};
    auto v2 = std::views::take(v1, 1);

    using iterator2_t = std::ranges::iterator_t<decltype(v2)>;
    EXPECT_TRUE(std::indirectly_readable<iterator2_t>);
    EXPECT_TRUE(ranges::indirectly_readable<iterator2_t>); // this failed
    EXPECT_TRUE(std::input_iterator<iterator2_t>);
    EXPECT_TRUE(ranges::input_iterator<iterator2_t>); // this failed
#endif                                                // __cpp_lib_ranges
}
