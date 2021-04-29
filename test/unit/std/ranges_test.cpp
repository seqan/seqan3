// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/ranges>
#include <string>

#include <range/v3/view/take.hpp>

TEST(ranges_test, combine_std_with_range_v3)
{
    std::string str{"foo"};
    auto take_first = str | std::views::take(5) | ranges::view::take(1);

    EXPECT_EQ(*std::ranges::begin(take_first), 'f');
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
#endif // __cpp_lib_ranges
}
