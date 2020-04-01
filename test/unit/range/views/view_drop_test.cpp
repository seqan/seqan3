// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <range/v3/view/unique.hpp>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor(3);
    EXPECT_EQ("bar", std::string(v));

    // function notation
    std::string v2{adaptor(vec, 3) | seqan3::views::to<std::string>};
    EXPECT_EQ("bar", v2);

    // combinability
    auto v3 = vec | adaptor(1) | adaptor(1) | ranges::view::unique;
    EXPECT_EQ("obar", v3 | seqan3::views::to<std::string>);
    std::string v3b = vec | std::views::reverse | adaptor(3) | ranges::view::unique | seqan3::views::to<std::string>;
    EXPECT_EQ("of", v3b);

    // store arg
    auto a0 = adaptor(3);
    auto v4 = vec | a0;
    EXPECT_EQ("bar", v4 | seqan3::views::to<std::string>);

    // store combined
    auto a1 = adaptor(1) | adaptor(1) | ranges::view::unique;
    auto v5 = vec | a1;
    EXPECT_EQ("obar", v5 | seqan3::views::to<std::string>);
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor)
{
    std::vector vec{1, 2, 3};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
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
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), int>));

    auto v2 = vec | seqan3::views::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v2), int>));
}

// ============================================================================
//  view_drop
// ============================================================================

TEST(view_drop, regular)
{
    do_test(seqan3::views::drop, "foobar");
}

TEST(view_drop, concepts)
{
    do_concepts(seqan3::views::drop(3));
}

TEST(view_drop, underlying_is_shorter)
{
    std::string vec{"foobar"};
    EXPECT_NO_THROW(( seqan3::views::drop(vec, 4) )); // no parsing

    std::string v;
    // full parsing on conversion
    EXPECT_NO_THROW(( v = vec
                        | seqan3::views::single_pass_input
                        | seqan3::views::drop(4)
                        | seqan3::views::to<std::string>));
    EXPECT_EQ("ar", v);
}

TEST(view_drop, type_erasure)
{
    {   // string overload
        std::string const urange{"foobar"};

        auto v = seqan3::views::drop(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange.substr(3,3))));
    }

    {   // stringview overload
        std::string_view urange{"foobar"};

        auto v = seqan3::views::drop(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange.substr(3,3))));
    }

    {   // contiguous overload
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::drop(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{4, 5, 6})));
    }

    {   // contiguous overload
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::drop(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{4, 5, 6})));
    }

    {   // generic overload (random-access container)
        std::deque<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::drop(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), std::ranges::subrange<typename std::deque<int>::iterator,
                                                                     typename std::deque<int>::iterator>>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{4, 5, 6})));
    }

    {   // no type erasure (bidirectional container)
        std::list<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3::views::drop(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), decltype(std::views::drop(urange, 3))>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{4, 5, 6})));
    }

    {   // no type erasure (input view)
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = urange | std::views::filter([] (int) { return true; });
        auto v2 = seqan3::views::drop(v, 3);

        EXPECT_TRUE((std::same_as<decltype(v2), decltype(std::views::drop(v, 3))>));
        EXPECT_TRUE((std::ranges::equal(v2, std::vector{4, 5, 6})));
    }
}
