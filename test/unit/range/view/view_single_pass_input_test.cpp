// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <gtest/gtest.h>

#include <vector>
#include <type_traits>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/istream_range.hpp>

#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/std/ranges>

template <typename rng_type>
class single_pass_input : public ::testing::Test
{

    virtual void SetUp()
    {
        data = get_data();
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
        else if constexpr (std::is_same_v<std::remove_cv_t<rng_type>, ranges::v3::istream_range<char>>)
        {
            return std::istringstream{"12345"};
        }
        else /*(std::is_same_v<rng_type, ranges::v3::istream_range<int>>)*/
        {
            return std::istringstream{"1 2 3 4 5"};
        }
    }

public:

    decltype(std::declval<single_pass_input>().get_data()) data;
};

// add all <out_rng,in_rng> pairs here.
using underlying_range_types = ::testing::Types<std::vector<char>,
                                                std::vector<int>,
                                                std::vector<char> const,
                                                ranges::istream_range<char>,
                                                ranges::istream_range<int>>;

TYPED_TEST_CASE(single_pass_input, underlying_range_types);

using namespace seqan3;

TYPED_TEST(single_pass_input, view_concept)
{
    using view_t = detail::single_pass_input_view<std::add_lvalue_reference_t<TypeParam>>;
    EXPECT_TRUE((std::is_base_of_v<ranges::view_base, view_t>));

    EXPECT_TRUE((std::Sentinel<sentinel_t<view_t>, iterator_t<view_t>>));
    EXPECT_TRUE(std::ranges::Range<view_t>);
    EXPECT_TRUE(std::ranges::View<view_t>);
    EXPECT_TRUE(std::ranges::InputRange<view_t>);
    EXPECT_FALSE(std::ranges::ForwardRange<view_t>);
}

TYPED_TEST(single_pass_input, view_construction)
{
    using view_t = detail::single_pass_input_view<std::add_lvalue_reference_t<TypeParam>>;
    EXPECT_TRUE(std::is_default_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<view_t>);
    EXPECT_TRUE(std::is_move_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<view_t>);
    EXPECT_TRUE(std::is_move_assignable_v<view_t>);
    EXPECT_TRUE(std::is_destructible_v<view_t>);

    {  // from lvalue
        TypeParam p{this->data};
        [[maybe_unused]] detail::single_pass_input_view v{p};
    }

    {  // from rvalue
        [[maybe_unused]] detail::single_pass_input_view v{TypeParam{this->data}};
    }
}

TYPED_TEST(single_pass_input, view_begin)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    using iterator_type = ranges::iterator_t<decltype(view)>;

    EXPECT_TRUE((std::is_same_v<decltype(view.begin()), iterator_type>));
    EXPECT_EQ(*view.begin(), *p.begin());
}

TYPED_TEST(single_pass_input, view_end)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    using sentinel_type = ranges::sentinel_t<decltype(view)>;

    EXPECT_TRUE((std::is_same_v<decltype(view.end()), sentinel_type>));
}

TYPED_TEST(single_pass_input, view_iterate)
{
    TypeParam p{this->data};

    detail::single_pass_input_view view{p};

    if constexpr (std::is_base_of_v<std::ios_base, decltype(this->data)>)
    {
        for (auto && elem : view)
        {
            EXPECT_EQ(elem, *p.begin());
        }
    }
    else
    {
        auto zipper = ranges::zip_view(p, view);
        for (auto it = zipper.begin(); it != zipper.end(); ++it)
        {
            EXPECT_EQ(std::get<0>(*it), std::get<1>(*it));
        }
    }
}

TYPED_TEST(single_pass_input, iterator_concepts)
{
    using view_type = detail::single_pass_input_view<std::add_lvalue_reference_t<TypeParam>>;
    EXPECT_TRUE((std::InputIterator<ranges::iterator_t<view_type>>));
    EXPECT_FALSE((std::ForwardIterator<ranges::iterator_t<view_type>>));
}

TYPED_TEST(single_pass_input, iterator_construction)
{
    using view_type = detail::single_pass_input_view<std::add_lvalue_reference_t<TypeParam>>;
    using iterator_type = ranges::iterator_t<view_type>;
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
    if constexpr (std::is_same_v<ranges::range_value_type_t<TypeParam>, char>)
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

    if constexpr (std::is_same_v<ranges::range_value_type_t<TypeParam>, char>)
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
    using view_type     = detail::single_pass_input_view<std::add_lvalue_reference_t<TypeParam>>;
    using iterator_type = iterator_t<view_type>;
    using sentinel_type = sentinel_t<view_type>;

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
    if constexpr (std::is_same_v<ranges::range_value_type_t<TypeParam>, char>)
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
    if constexpr (std::is_same_v<ranges::range_value_type_t<TypeParam>, char>)
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
