// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

TEST(view_persist, delegate_to_view_all)
{
    std::string vec{"foo"};

    // pipe notation
    auto v = vec | view::persist;
    EXPECT_EQ("foo", std::string{v});

    // function notation
    std::string v2 = view::persist(vec);
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | view::persist | ranges::view::unique;
    EXPECT_EQ("fo", std::string{v3});
    std::string v3b = vec | std::view::reverse | view::persist | ranges::view::unique;
    EXPECT_EQ("of", v3b);
}

TEST(view_persist, wrap_temporary)
{
    // pipe notation
    auto v = std::string{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string(v));

    // function notation
    std::string v2 = view::persist(std::string{"foo"});
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = std::string{"foo"} | view::persist | ranges::view::unique;
    EXPECT_EQ("fo", std::string(v3));
    std::string v3b = std::string{"foo"} | view::persist | std::view::filter(is_char<'o'>) | ranges::view::unique;
    EXPECT_EQ("o", v3b);
}

TEST(view_persist, const_)
{
    // inner const
    using t = std::string const;
    auto v = t{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string{v});

    // outer const
    auto const & v2 = std::string{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string{v2});

    // inner + outer const
    using t = std::string const;
    auto const & v3 = t{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string{v3});
}

TEST(view_persist, concepts)
{
    using namespace seqan3::test;
    auto v1 = std::string{"foo"} | view::persist;

    EXPECT_TRUE((preserved<std::string, decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Contiguous, Sized, 
                                                         Common, Output, ConstIterable})));
    EXPECT_TRUE((guaranteed<std::string, decltype(v1)>({View, Viewable})));
    EXPECT_TRUE((weak_guaranteed<decltype(v1)>({Viewable})));
}
