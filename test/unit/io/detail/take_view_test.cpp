// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <deque>
#include <list>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/io/detail/take_view.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

inline auto constexpr seqan3_views_take = seqan3::detail::take_fn<false, false>{};

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    using namespace std::literals;

    // pipe notation
    EXPECT_RANGE_EQ("foo"sv, vec | adaptor(3) );

    // function notation
    EXPECT_RANGE_EQ("foo"sv, adaptor(vec, 3));

    // combinability
    EXPECT_RANGE_EQ("fo"sv, vec | adaptor(3) | adaptor(3) | std::views::take(2));
    EXPECT_RANGE_EQ("rab"sv, vec | std::views::reverse | adaptor(3) | std::views::take(3));
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const exactly)
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

    auto v3 = vec | std::views::transform([] (auto && v) { return v; }) | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::view<decltype(v3)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v3)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v3)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v3), int>));

    auto v2 = vec | seqan3::views::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(v2)>, exactly);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v2), int>));

    // explicit test for non const-iterable views
    // https://github.com/seqan/seqan3/pull/1734#discussion_r408829267
    auto const & v2_cref = v2;

    EXPECT_FALSE(std::ranges::input_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::view<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2_cref)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2_cref)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v2_cref)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2_cref), int>));
}

// ============================================================================
//  view_take
// ============================================================================

TEST(view_take, regular)
{
    do_test(seqan3_views_take, "foobar");
}

TEST(view_take, concepts)
{
    do_concepts(seqan3_views_take(3), false);
}

TEST(view_take, underlying_is_shorter)
{
    using namespace std::literals;

    std::string vec{"foo"};
    EXPECT_NO_THROW(( seqan3_views_take(vec, 4) )); // no parsing

    std::string v;
    // full parsing on conversion
    EXPECT_RANGE_EQ("foo"sv, vec | seqan3::views::single_pass_input | seqan3_views_take(4));
}

TEST(view_take, type_erasure)
{
    {   // string const overload
        std::string const urange{"foobar"};

        auto v = seqan3_views_take(urange, 3);

        EXPECT_SAME_TYPE(decltype(v), std::string_view);
        EXPECT_RANGE_EQ(v, urange.substr(0,3));
    }

    {   // stringview overload
        std::string_view urange{"foobar"};

        auto v = seqan3_views_take(urange, 3);

        EXPECT_SAME_TYPE(decltype(v), std::string_view);
        EXPECT_RANGE_EQ(v, urange.substr(0,3));
    }

    {   // contiguous overload
        std::vector<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3_views_take(urange, 3);

        EXPECT_SAME_TYPE(decltype(v), (std::span<int, std::dynamic_extent>));
        EXPECT_RANGE_EQ(v, (std::vector{1, 2, 3}));
    }

    {   // contiguous overload
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3_views_take(urange, 3);

        EXPECT_SAME_TYPE(decltype(v), (std::span<int, std::dynamic_extent>));
        EXPECT_RANGE_EQ(v, (std::vector{1, 2, 3}));
    }

    {   // random-access overload
        std::deque<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3_views_take(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), std::ranges::subrange<typename std::deque<int>::iterator,
                                                                     typename std::deque<int>::iterator>>));
        EXPECT_RANGE_EQ(v, (std::vector{1, 2, 3}));
    }

    {   // generic overload (bidirectional container)
        std::list<int> urange{1, 2, 3, 4, 5, 6};

        auto v = seqan3_views_take(urange, 3);

        EXPECT_TRUE((std::same_as<decltype(v), seqan3::detail::view_take<std::views::all_t<std::list<int> &>,
                                                                         false, false>>));
        EXPECT_RANGE_EQ(v, (std::vector{1, 2, 3}));
    }

    {   // generic overload (view)
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v = urange | std::views::filter([] (int) { return true; });
        auto v2 = seqan3_views_take(v, 3);

        EXPECT_SAME_TYPE(decltype(v2), (seqan3::detail::view_take<decltype(v), false, false>));
        EXPECT_RANGE_EQ(v2, (std::vector{1, 2, 3}));
    }

    {   // generic overload (random access, non-sized, pointer as iterator)
        std::array<int, 6> urange{1, 2, 3, 4, 5, 6};

        auto v0 = std::span{urange};
        auto v1 = v0 | std::views::take_while([] (int i) { return i < 6; });
        auto v2 = seqan3_views_take(v1, 3);

        EXPECT_SAME_TYPE(decltype(v2), (seqan3::detail::view_take<decltype(v1), false, false>));
        EXPECT_RANGE_EQ(v2, (std::vector{1, 2, 3}));
    }
}
