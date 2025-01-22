// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <concepts>
#include <deque>
#include <iostream>
#include <list>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

// ============================================================================
//  test templates
// ============================================================================

TEST(type_reduce, string_overload)
{
    {
        std::string urange{"foobar"};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_FALSE((std::same_as<decltype(v), std::string_view>)); // only returns string_view for string const
        EXPECT_RANGE_EQ(v, urange);
    }

    {
        std::string urange_{"foobar"};
        std::string_view urange{urange_};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_SAME_TYPE(decltype(v), std::string_view);
        EXPECT_RANGE_EQ(v, urange);
    }

    {
        std::string const urange{"foobar"};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_SAME_TYPE(decltype(v), std::string_view);
        EXPECT_RANGE_EQ(v, urange);
    }
}

TEST(type_reduce, contiguous_overload)
{
    {
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_SAME_TYPE(decltype(v), (std::span<int, std::dynamic_extent>));
        EXPECT_RANGE_EQ(v, urange);
    }

    {
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_SAME_TYPE(decltype(v), (std::span<int, std::dynamic_extent>));
        EXPECT_RANGE_EQ(v, urange);
    }
}

TEST(type_reduce, random_access_overload)
{
    std::deque<int> urange{1, 2, 3, 4, 5, 6};

    auto v = seqan3::views::type_reduce(urange);

    EXPECT_TRUE(
        (std::same_as<decltype(v),
                      std::ranges::subrange<typename std::deque<int>::iterator, typename std::deque<int>::iterator>>));
    EXPECT_RANGE_EQ(v, urange);
}

TEST(type_reduce, generic_overload)
{
    { // bidirectional container
        std::list<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::type_reduce(urange);

        EXPECT_SAME_TYPE(decltype(v), std::views::all_t<std::list<int> &>);
        EXPECT_RANGE_EQ(v, urange);
    }

    { // view
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = urange
               | std::views::filter(
                     [](int)
                     {
                         return true;
                     });
        auto v2 = seqan3::views::type_reduce(v);

        EXPECT_SAME_TYPE(decltype(v2), std::views::all_t<decltype(v)>);
        EXPECT_RANGE_EQ(v2, urange);
    }
}
