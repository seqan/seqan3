// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>
#include <deque>
#include <list>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

TEST(view_all, string_overload)
{
    {
        std::string urange{"foobar"};

        auto v = view::all(urange);

        EXPECT_FALSE((std::Same<decltype(v), std::string_view>)); // only returns string_view for string const
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {
        std::string urange_{"foobar"};
        std::string_view urange{urange_};

        auto v = view::all(urange);

        EXPECT_TRUE((std::Same<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {
        std::string const urange{"foobar"};

        auto v = view::all(urange);

        EXPECT_TRUE((std::Same<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }
}

TEST(view_all, contiguous_overload)
{
    {
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = view::all(urange);

        EXPECT_TRUE((std::Same<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = view::all(urange);

        EXPECT_TRUE((std::Same<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }
}

TEST(view_all, random_access_overload)
{
    std::deque<int> urange{1, 2, 3, 4, 5, 6};

    auto v = view::all(urange);

    EXPECT_TRUE((std::Same<decltype(v),
                          std::ranges::subrange<typename std::deque<int>::iterator, typename std::deque<int>::iterator>>));
    EXPECT_TRUE((std::ranges::equal(v, urange)));
}

TEST(view_all, generic_overload)
{
    {   // bidirectional container
        std::list<int> urange{1, 2, 3, 4, 5, 6};

        auto v = view::all(urange);

        EXPECT_TRUE((std::Same<decltype(v), std::ranges::all_view<std::list<int> &>>));
        EXPECT_TRUE((std::ranges::equal(v, urange)));
    }

    {   // view
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = urange | std::view::filter([] (int) { return true; });
        auto v2 = view::all(v);

        EXPECT_TRUE((std::Same<decltype(v2), std::ranges::all_view<decltype(v)>>));
        EXPECT_TRUE((std::ranges::equal(v2, urange)));
    }
}
