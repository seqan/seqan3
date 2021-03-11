// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <string>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/test/expect_range_eq.hpp>

// ============================================================================
//  test templates
// ============================================================================

TEST(view_persist, delegate_to_view_all)
{
    using namespace std::literals;

    std::string vec{"foo"};

    // pipe notation
    auto v = vec | seqan3::views::persist;
    EXPECT_RANGE_EQ("foo"sv, v);

    // function notation
    EXPECT_RANGE_EQ("foo"sv, seqan3::views::persist(vec));

    // combinability
    EXPECT_RANGE_EQ("fo"sv, vec | seqan3::views::persist | ranges::views::unique);
    EXPECT_RANGE_EQ("of"sv, vec | std::views::reverse | seqan3::views::persist | ranges::views::unique);

    // store combined
    auto a1 = seqan3::views::persist | ranges::views::unique;
    EXPECT_RANGE_EQ("fo"sv, vec | a1);
}

TEST(view_persist, wrap_temporary)
{
    using namespace std::literals;

    // pipe notation
    EXPECT_RANGE_EQ("foo"sv, std::string{"foo"} | seqan3::views::persist);

    // function notation
    EXPECT_RANGE_EQ("foo"sv, seqan3::views::persist(std::string{"foo"}));

    // combinability
    EXPECT_RANGE_EQ("fo"sv, std::string{"foo"} | seqan3::views::persist | ranges::views::unique);
    EXPECT_RANGE_EQ("o"sv, std::string{"foo"} | seqan3::views::persist
                                              | std::views::filter([](char const chr){return chr == 'o';})
                                              | ranges::views::unique);
}

TEST(view_persist, const)
{
    using namespace std::literals;

    // inner const
    using t = std::string const;
    EXPECT_RANGE_EQ("foo"sv, t{"foo"} | seqan3::views::persist);

    // outer const
    auto const & v2 = std::string{"foo"} | seqan3::views::persist;
    EXPECT_RANGE_EQ("foo"sv, v2);

    // inner + outer const
    using t = std::string const;
    auto const & v3 = t{"foo"} | seqan3::views::persist;
    EXPECT_RANGE_EQ("foo"sv, v3);
}

TEST(view_persist, concepts)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(std::string{"foo"})>);
    EXPECT_FALSE(std::ranges::view<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::common_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(std::string{"foo"})>);
    EXPECT_TRUE((std::ranges::output_range<decltype(std::string{"foo"}), char>));

    auto v1 = std::string{"foo"} | seqan3::views::persist;

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));
}
