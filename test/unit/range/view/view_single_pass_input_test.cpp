// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>
#include <type_traits>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/istream_range.hpp>

#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/std/ranges>

template <typename rng_type>
class single_pass_input : public ::testing::Test
{

    virtual void SetUp()
    {
        data = get_data();
        cmp_data = get_data();
    }

    auto get_data()
    {
        if constexpr (std::is_same_v<std::remove_cv_t<rng_type>, std::vector<char>>)
        {
            return rng_type{'1','2','3','4','5'};
        }
        else if constexpr (std::is_same_v<std::remove_cv_t<rng_type>, std::vector<int>>)
        {
            return rng_type{1, 2, 3, 4, 5};
        }
        else if constexpr (std::is_same_v<std::remove_cv_t<rng_type>, ranges::istream_view<char>>)
        {
            return std::istringstream{"12345"};
        }
        else if constexpr (std::is_same_v<rng_type, ranges::istream_view<int>>)
        {
            return std::istringstream{"1 2 3 4 5"};
        }
    }

public:

    decltype(std::declval<single_pass_input>().get_data()) data;
    decltype(std::declval<single_pass_input>().get_data()) cmp_data;
};

// add all <out_rng,in_rng> pairs here.
using underlying_range_types = ::testing::Types<std::vector<char>,
                                                std::vector<int>,
                                                std::vector<char> const,
                                                ranges::istream_view<char>,
                                                ranges::istream_view<int>>;

TYPED_TEST_CASE(single_pass_input, underlying_range_types);

using namespace seqan3;

TYPED_TEST(single_pass_input, view_concept)
{
    using rng_t = decltype(std::declval<TypeParam>() | view::persist);
    using view_t = detail::single_pass_input_view<rng_t>;
    EXPECT_TRUE((std::is_base_of_v<ranges::view_base, view_t>));
    EXPECT_TRUE((std::Sentinel<std::ranges::sentinel_t<view_t>, std::ranges::iterator_t<view_t>>));
    EXPECT_TRUE(std::ranges::Range<view_t>);
    EXPECT_TRUE(std::ranges::View<view_t>);
    EXPECT_TRUE(std::ranges::InputRange<view_t>);
    EXPECT_FALSE(std::ranges::ForwardRange<view_t>);
}

TYPED_TEST(single_pass_input, view_construction)
{
    using rng_t = decltype(std::declval<TypeParam>() | view::persist);
    using view_t = detail::single_pass_input_view<rng_t>;
    EXPECT_TRUE(std::is_default_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<view_t>);
    EXPECT_TRUE(std::is_move_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<view_t>);
    EXPECT_TRUE(std::is_move_assignable_v<view_t>);
    EXPECT_TRUE(std::is_destructible_v<view_t>);

    {  // from lvalue container
        TypeParam p{this->data};
        [[maybe_unused]] detail::single_pass_input_view v{p};
        EXPECT_TRUE((std::is_same_v<decltype(v), detail::single_pass_input_view<std::remove_reference_t<decltype(p | view::all)>>>));
    }

    {  // from view
        [[maybe_unused]] detail::single_pass_input_view v{TypeParam{this->data} | view::persist};
        EXPECT_TRUE((std::is_same_v<decltype(v),
                    detail::single_pass_input_view<decltype(TypeParam{this->data} | view::persist)>>));
    }
}

TYPED_TEST(single_pass_input, view_begin)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    using iterator_type = std::ranges::iterator_t<decltype(view)>;

    EXPECT_TRUE((std::is_same_v<decltype(view.begin()), iterator_type>));
    EXPECT_EQ(*view.begin(), *p.begin());
}

TYPED_TEST(single_pass_input, view_end)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    using sentinel_type = std::ranges::sentinel_t<decltype(view)>;

    EXPECT_TRUE((std::is_same_v<decltype(view.end()), sentinel_type>));
}

TYPED_TEST(single_pass_input, view_iterate)
{
    TypeParam p{this->data};

    if constexpr (std::is_base_of_v<std::ios_base, decltype(this->data)>)
    {
        // Single pass input is only movable.
        detail::single_pass_input_view view{std::move(p)};

        TypeParam tmp{this->cmp_data};
        auto tmp_it = tmp.begin();
        for (auto && elem : view)
        {
            EXPECT_EQ(elem, *tmp_it);
            ++tmp_it;
        }
    }
    else
    {
        detail::single_pass_input_view view{p};
        TypeParam tmp{this->cmp_data};
        auto zipper = ranges::view::zip(tmp, std::move(view));
        for (auto it = zipper.begin(); it != zipper.end(); ++it)
        {
            EXPECT_EQ(std::get<0>(*it), std::get<1>(*it));
        }
    }
}

TYPED_TEST(single_pass_input, iterator_concepts)
{
    using view_type = detail::single_pass_input_view<decltype(std::declval<TypeParam>() | view::persist)>;
    EXPECT_TRUE((std::InputIterator<std::ranges::iterator_t<view_type>>));
    EXPECT_FALSE((std::ForwardIterator<std::ranges::iterator_t<view_type>>));
}

TYPED_TEST(single_pass_input, iterator_construction)
{
    using view_type = detail::single_pass_input_view<decltype(std::declval<TypeParam>() | view::persist)>;
    using iterator_type = std::ranges::iterator_t<view_type>;
    EXPECT_TRUE(std::is_default_constructible_v<iterator_type>);
    EXPECT_TRUE(std::is_copy_constructible_v<iterator_type>);
    EXPECT_TRUE(std::is_move_constructible_v<iterator_type>);
    EXPECT_TRUE(std::is_copy_assignable_v<iterator_type>);
    EXPECT_TRUE(std::is_move_assignable_v<iterator_type>);
    EXPECT_TRUE(std::is_destructible_v<iterator_type>);
}

TYPED_TEST(single_pass_input, iterator_pre_increment)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    auto it = view.begin();
    if constexpr (std::is_same_v<value_type_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it,   '1');
        EXPECT_EQ(*++it, '2');
        EXPECT_EQ(*++it, '3');
        EXPECT_EQ(*++it, '4');
        EXPECT_EQ(*++it, '5');
    }
    else
    {
        EXPECT_EQ(*it,   1);
        EXPECT_EQ(*++it, 2);
        EXPECT_EQ(*++it, 3);
        EXPECT_EQ(*++it, 4);
        EXPECT_EQ(*++it, 5);
    }
}

TYPED_TEST(single_pass_input, iterator_post_increment)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    auto it = view.begin();
    EXPECT_TRUE((std::Same<decltype(it++), void>));

    if constexpr (std::is_same_v<value_type_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it, '1');
        it++;
        EXPECT_EQ(*it, '2');
        it++;
        EXPECT_EQ(*it, '3');
        it++;
        EXPECT_EQ(*it, '4');
        it++;
        EXPECT_EQ(*it, '5');
    }
    else
    {
        EXPECT_EQ(*it, 1);
        it++;
        EXPECT_EQ(*it, 2);
        it++;
        EXPECT_EQ(*it, 3);
        it++;
        EXPECT_EQ(*it, 4);
        it++;
        EXPECT_EQ(*it, 5);
    }
}

TYPED_TEST(single_pass_input, iterator_eq_comparison)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};
    EXPECT_FALSE(view.begin() == view.end());

    auto it = view.begin();
    ++it;
    ++it;
    ++it;
    ++it;
    EXPECT_FALSE(view.begin() == view.end());
    ++it;
    EXPECT_TRUE(view.begin() == view.end());
}

TYPED_TEST(single_pass_input, iterator_neq_comparison)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};
    EXPECT_TRUE(view.begin() != view.end());

    auto it = view.begin();
    ++it;
    ++it;
    ++it;
    ++it;
    EXPECT_TRUE(view.begin() != view.end());
    ++it;
    EXPECT_FALSE(view.begin() != view.end());
}

TYPED_TEST(single_pass_input, sentinel_concepts)
{
    using view_type     = detail::single_pass_input_view<decltype(std::declval<TypeParam>() | view::persist)>;
    using iterator_type = std::ranges::iterator_t<view_type>;
    using sentinel_type = std::ranges::sentinel_t<view_type>;

    EXPECT_TRUE((std::Sentinel<sentinel_type, iterator_type>));
    EXPECT_FALSE((std::SizedSentinel<sentinel_type, iterator_type>));
}

TYPED_TEST(single_pass_input, sentinel_eq_comparison)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};
    EXPECT_FALSE(view.end() == view.begin());

    auto it = view.begin();
    ++it;
    ++it;
    ++it;
    ++it;
    EXPECT_FALSE(view.end() == view.begin());
    ++it;
    EXPECT_TRUE(view.end() == view.begin());
}

TYPED_TEST(single_pass_input, sentinel_neq_comparison)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};
    EXPECT_TRUE(view.end() != view.begin());

    auto it = view.begin();
    ++it;
    ++it;
    ++it;
    ++it;
    EXPECT_TRUE(view.end() != view.begin());
    ++it;
    EXPECT_FALSE(view.end() != view.begin());
}

TYPED_TEST(single_pass_input, fn_functional)
{
    // use case 1: functional;
    TypeParam p{this->data};

    auto view = ranges::view::take(view::single_pass_input(p), 3);

    auto it = view.begin();
    if constexpr (std::is_same_v<value_type_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it,   '1');
        EXPECT_EQ(*++it, '2');
        EXPECT_EQ(*++it, '3');
    }
    else
    {
        EXPECT_EQ(*it,   1);
        EXPECT_EQ(*++it, 2);
        EXPECT_EQ(*++it, 3);
    }
    ++it;
    EXPECT_TRUE(view.end() == it);
}

TYPED_TEST(single_pass_input, fn_pipeable)
{
    TypeParam p{this->data};

    auto view = p | view::single_pass_input | ranges::view::take(3);
    auto it = view.begin();
    if constexpr (std::is_same_v<value_type_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it,   '1');
        EXPECT_EQ(*++it, '2');
        EXPECT_EQ(*++it, '3');
    }
    else
    {
        EXPECT_EQ(*it,   1);
        EXPECT_EQ(*++it, 2);
        EXPECT_EQ(*++it, 3);
    }
    ++it;
    EXPECT_TRUE(view.end() == it);
}
