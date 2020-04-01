// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>
#include <deque>
#include <list>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

// ============================================================================
//  test templates
// ============================================================================

TEST(type_reduce, string_overload)
{
    {
        std::string urange{"foobar"};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_FALSE((std::same_as<decltype(v), std::string_view>)); // only returns string_view for string const
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {
        std::string urange_{"foobar"};
        std::string_view urange{urange_};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_TRUE((std::same_as<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {
        std::string const urange{"foobar"};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_TRUE((std::same_as<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }
}

TEST(type_reduce, contiguous_overload)
{
    {
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_TRUE((std::same_as<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_TRUE((std::same_as<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }
}

TEST(type_reduce, random_access_overload)
{
    std::deque<int> urange{1, 2, 3, 4, 5, 6};

    auto v = seqan3::views::type_reduce(urange);

    EXPECT_TRUE((std::same_as<decltype(v), std::ranges::subrange<typename std::deque<int>::iterator,
                                                                 typename std::deque<int>::iterator>>));
    EXPECT_TRUE((std::ranges::equal(v, urange)));
}

TEST(type_reduce, generic_overload)
{
    {   // bidirectional container
        std::list<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_TRUE((std::same_as<decltype(v), std::ranges::all_view<std::list<int> &>>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {   // view
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = urange | std::views::filter([] (int) { return true; });
        auto v2 = seqan3::views::type_reduce(v);

        EXPECT_TRUE((std::same_as<decltype(v2), std::ranges::all_view<decltype(v)>>));
        EXPECT_TRUE((std::ranges::equal(v2, urange)));
    }
}
