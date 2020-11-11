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

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

// ============================================================================
//  test templates
// ============================================================================

TEST(view_persist, delegate_to_view_all)
{
    std::string vec{"foo"};

    // pipe notation
    auto v = vec | seqan3::views::persist;
    EXPECT_EQ("foo", v | seqan3::views::to<std::string>);

    // function notation
    std::string v2 = seqan3::views::persist(vec) | seqan3::views::to<std::string>;
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | seqan3::views::persist | ranges::views::unique;
    EXPECT_EQ("fo", v3 | seqan3::views::to<std::string>);
    std::string v3b = vec
                    | std::views::reverse
                    | seqan3::views::persist
                    | ranges::views::unique
                    | seqan3::views::to<std::string>;
    EXPECT_EQ("of", v3b);

    // store combined
    auto a1 = seqan3::views::persist | ranges::views::unique;
    auto v5 = vec | a1;
    EXPECT_EQ("fo", v5 | seqan3::views::to<std::string>);
}

TEST(view_persist, wrap_temporary)
{
    // pipe notation
    auto v = std::string{"foo"} | seqan3::views::persist;
    EXPECT_EQ("foo", v | seqan3::views::to<std::string>);

    // function notation
    std::string v2 = seqan3::views::persist(std::string{"foo"}) | seqan3::views::to<std::string>;
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = std::string{"foo"} | seqan3::views::persist | ranges::views::unique;
    EXPECT_EQ("fo", v3 | seqan3::views::to<std::string>);
    std::string v3b = std::string{"foo"}
                    | seqan3::views::persist
                    | std::views::filter(seqan3::is_char<'o'>)
                    | ranges::views::unique
                    | seqan3::views::to<std::string>;
    EXPECT_EQ("o", v3b);
}

TEST(view_persist, const)
{
    // inner const
    using t = std::string const;
    auto v = t{"foo"} | seqan3::views::persist;
    EXPECT_EQ("foo", v | seqan3::views::to<std::string>);

    // outer const
    auto const & v2 = std::string{"foo"} | seqan3::views::persist;
    EXPECT_EQ("foo", v2 | seqan3::views::to<std::string>);

    // inner + outer const
    using t = std::string const;
    auto const & v3 = t{"foo"} | seqan3::views::persist;
    EXPECT_EQ("foo", v3 | seqan3::views::to<std::string>);
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
