// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

#include <seqan3/range/view/slice.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor(1, 4);
    EXPECT_EQ("oob", std::string(v));

    // function notation
    std::string v2{adaptor(vec, 1, 4)};
    EXPECT_EQ("oob", v2);

    // combinability
    auto v3 = vec | adaptor(0, 4) | adaptor(1, 3) | ranges::view::unique;
    EXPECT_EQ("o", std::string(v3));
    std::string v3b = vec | std::view::reverse | adaptor(1, 4) | ranges::view::unique;
    EXPECT_EQ("abo", v3b);

    // store arg
    auto a0 = adaptor(1, 4);
    auto v4 = vec | a0;
    EXPECT_EQ("oob", std::string(v4));

    // store combined
    auto a1 = adaptor(0, 4) | adaptor(1, 3) | ranges::view::unique;
    auto v5 = vec | a1;
    EXPECT_EQ("o", std::string(v5));

    // store combined in middle
    auto a2 = std::view::reverse | adaptor(1, 4) | ranges::view::unique;
    auto v6 = vec | a2;
    EXPECT_EQ("abo", std::string(v6));
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const exactly)
{
    std::vector vec{1, 2, 3};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(ConstIterableRange<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), int>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(ConstIterableRange<decltype(v1)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), int>));

    auto v2 = vec | view::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::ForwardRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::BidirectionalRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::RandomAccessRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::View<decltype(v2)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(v2)>, exactly);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v2)>);
    EXPECT_FALSE(ConstIterableRange<decltype(v2)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v2), int>));
}

// ============================================================================
//  view_slice
// ============================================================================

TEST(view_slice, regular)
{
    do_test(view::slice, "foobar");
}

TEST(view_slice, concepts)
{
    do_concepts(view::slice(1, 4), false);
}

TEST(view_slice, underlying_is_shorter)
{
    std::string vec{"foobar"};
    EXPECT_NO_THROW(( view::slice(vec, 1, 4) )); // no parsing

    std::string v;
    EXPECT_NO_THROW(( v = vec | view::single_pass_input | view::slice(1, 4) )); // full parsing on conversion
    EXPECT_EQ("oob", v);
}

TEST(view_slice, type_erasure)
{
    {   // string overload
        std::string const urange{"foobar"};

        auto v = view::slice(urange, 1, 4);

        EXPECT_TRUE((std::Same<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange.substr(1,3))));
    }

    {   // stringview overload
        std::string_view urange{"foobar"};

        auto v = view::slice(urange, 1, 4);

        EXPECT_TRUE((std::Same<decltype(v), std::string_view>));
        EXPECT_TRUE((std::ranges::equal(v, urange.substr(1,3))));
    }

    {   // contiguous overload
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = view::slice(urange, 1, 4);

        EXPECT_TRUE((std::Same<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{2, 3, 4})));
    }

    {   // contiguous overload
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = view::slice(urange, 1, 4);

        EXPECT_TRUE((std::Same<decltype(v), std::span<int, std::dynamic_extent>>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{2, 3, 4})));
    }

    {   // random-access overload
        std::deque<int> urange{1, 2, 3, 4, 5, 6};

        auto v = view::slice(urange, 1, 4);

        EXPECT_TRUE((std::Same<decltype(v), std::ranges::subrange<typename std::deque<int>::iterator,
                                                                  typename std::deque<int>::iterator>>));
        EXPECT_TRUE((std::ranges::equal(v, std::vector{2, 3, 4})));
    }
}
