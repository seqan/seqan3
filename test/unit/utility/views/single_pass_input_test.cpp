// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ranges>
#include <type_traits>
#include <vector>

#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

template <typename rng_type>
class single_pass_input : public ::testing::Test
{

    virtual void SetUp()
    {
        data = get_data();
        cmp_data = get_data();
    }

    static auto get_data()
    {
        if constexpr (std::is_same_v<std::remove_cv_t<rng_type>, std::vector<char>>)
        {
            return rng_type{'1', '2', '3', '4', '5'};
        }
        else if constexpr (std::is_same_v<std::remove_cv_t<rng_type>, std::vector<int>>)
        {
            return rng_type{1, 2, 3, 4, 5};
        }
        else if constexpr (std::is_same_v<std::remove_cv_t<rng_type>,
                                          std::ranges::basic_istream_view<char, char, std::char_traits<char>>>)
        {
            return std::istringstream{"12345"};
        }
        else if constexpr (std::is_same_v<rng_type, std::ranges::basic_istream_view<int, char, std::char_traits<char>>>)
        {
            return std::istringstream{"1 2 3 4 5"};
        }
    }

public:
    decltype(get_data()) data;
    decltype(get_data()) cmp_data;
};

// add all <out_rng,in_rng> pairs here.
using underlying_range_types = ::testing::Types<std::vector<char>,
                                                std::vector<int>,
                                                std::ranges::basic_istream_view<char, char, std::char_traits<char>>,
                                                std::ranges::basic_istream_view<int, char, std::char_traits<char>>>;

TYPED_TEST_SUITE(single_pass_input, underlying_range_types, );

TYPED_TEST(single_pass_input, view_concept)
{
    using view_t = seqan3::detail::single_pass_input_view<std::views::all_t<TypeParam &>>;
    EXPECT_TRUE((std::derived_from<view_t, std::ranges::view_interface<view_t>>));
    EXPECT_TRUE((std::sentinel_for<std::ranges::sentinel_t<view_t>, std::ranges::iterator_t<view_t>>));
    EXPECT_TRUE(std::ranges::range<view_t>);
    EXPECT_TRUE(std::ranges::view<view_t>);
    EXPECT_TRUE(std::ranges::input_range<view_t>);
    EXPECT_FALSE((std::ranges::output_range<view_t, std::ranges::range_reference_t<view_t>>));
    EXPECT_FALSE(std::ranges::common_range<view_t>);
    EXPECT_FALSE(std::ranges::forward_range<view_t>);
    EXPECT_FALSE(std::ranges::bidirectional_range<view_t>);
    EXPECT_FALSE(std::ranges::random_access_range<view_t>);
}

TYPED_TEST(single_pass_input, deduction_guide_lvalue)
{
    TypeParam data_container{this->data};

    EXPECT_TRUE((std::ranges::viewable_range<TypeParam &>));

    // this is
    std::views::all_t<TypeParam &> data_view{data_container};
    seqan3::detail::single_pass_input_view<decltype(data_view)> v1{data_view};

    // same as
    seqan3::detail::single_pass_input_view v2{data_container};

    EXPECT_SAME_TYPE(decltype(v1), decltype(v2));
}

TYPED_TEST(single_pass_input, deduction_guide_view)
{
    TypeParam data_container{this->data};
    auto data_view = std::views::transform(data_container,
                                           [](auto const & in)
                                           {
                                               return in;
                                           });

    using uview_t = decltype(data_view);
    EXPECT_TRUE((std::ranges::viewable_range<uview_t>));
    EXPECT_TRUE((std::ranges::view<uview_t>));

    // this is
    std::views::all_t<uview_t> data_all_view{data_view};

    // same as
    seqan3::detail::single_pass_input_view<decltype(data_all_view)> v1{data_all_view};
    seqan3::detail::single_pass_input_view v2{std::move(data_view)};

    EXPECT_SAME_TYPE(decltype(v1), decltype(v2));
}

TYPED_TEST(single_pass_input, view_construction)
{
    using view_t = seqan3::detail::single_pass_input_view<std::views::all_t<TypeParam>>;
    EXPECT_TRUE(std::is_default_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<view_t>);
    EXPECT_TRUE(std::is_move_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<view_t>);
    EXPECT_TRUE(std::is_move_assignable_v<view_t>);
    EXPECT_TRUE(std::is_destructible_v<view_t>);

    { // from lvalue container
        TypeParam p{this->data};
        [[maybe_unused]] seqan3::detail::single_pass_input_view v{p};
    }

    { // from view
        [[maybe_unused]] seqan3::detail::single_pass_input_view v{TypeParam{this->data} | std::views::all};
    }
}

TYPED_TEST(single_pass_input, view_begin)
{
    using value_t = std::ranges::range_value_t<TypeParam>;
    TypeParam p{this->data};

    seqan3::detail::single_pass_input_view view{p};

    using iterator_type = std::ranges::iterator_t<decltype(view)>;

    EXPECT_TRUE((std::is_same_v<decltype(view.begin()), iterator_type>));

    value_t first_value = std::is_same_v<value_t, char> ? '1' : 1;
    EXPECT_EQ(*view.begin(), first_value);
}

TYPED_TEST(single_pass_input, view_end)
{
    TypeParam p{this->data};

    seqan3::detail::single_pass_input_view view{p};

    using sentinel_type = std::ranges::sentinel_t<decltype(view)>;

    EXPECT_TRUE((std::is_same_v<decltype(view.end()), sentinel_type>));
}

TYPED_TEST(single_pass_input, view_iterate)
{
    if constexpr (std::is_base_of_v<std::ios_base, decltype(this->data)>)
    {
        TypeParam p{this->data};

        // Single pass input is only movable.
        seqan3::detail::single_pass_input_view view{std::move(p)};

        TypeParam tmp{this->cmp_data};
        auto tmp_it = tmp.begin();
        for (auto && elem : view)
        {
            EXPECT_EQ(elem, *tmp_it);
            ++tmp_it;
        }
    }
}

TYPED_TEST(single_pass_input, iterator_concepts)
{
    using view_type = seqan3::detail::single_pass_input_view<std::views::all_t<TypeParam>>;
    EXPECT_TRUE((std::input_iterator<std::ranges::iterator_t<view_type>>));
    EXPECT_FALSE((std::forward_iterator<std::ranges::iterator_t<view_type>>));
}

TYPED_TEST(single_pass_input, iterator_construction)
{
    using view_type = seqan3::detail::single_pass_input_view<std::views::all_t<TypeParam>>;
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

    seqan3::detail::single_pass_input_view view{p};

    auto it = view.begin();
    if constexpr (std::is_same_v<std::ranges::range_value_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it, '1');
        EXPECT_EQ(*++it, '2');
        EXPECT_EQ(*++it, '3');
        EXPECT_EQ(*++it, '4');
        EXPECT_EQ(*++it, '5');
    }
    else
    {
        EXPECT_EQ(*it, 1);
        EXPECT_EQ(*++it, 2);
        EXPECT_EQ(*++it, 3);
        EXPECT_EQ(*++it, 4);
        EXPECT_EQ(*++it, 5);
    }
}

TYPED_TEST(single_pass_input, iterator_post_increment)
{
    TypeParam p{this->data};

    seqan3::detail::single_pass_input_view view{p};

    auto it = view.begin();

    if constexpr (std::is_same_v<std::ranges::range_value_t<TypeParam>, char>)
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

    seqan3::detail::single_pass_input_view view{p};
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

    seqan3::detail::single_pass_input_view view{p};
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
    using view_type = seqan3::detail::single_pass_input_view<std::views::all_t<TypeParam>>;
    using iterator_type = std::ranges::iterator_t<view_type>;
    using sentinel_type = std::ranges::sentinel_t<view_type>;

    EXPECT_TRUE((std::sentinel_for<sentinel_type, iterator_type>));
    EXPECT_FALSE((std::sized_sentinel_for<sentinel_type, iterator_type>));
}

TYPED_TEST(single_pass_input, sentinel_eq_comparison)
{
    TypeParam p{this->data};

    seqan3::detail::single_pass_input_view view{p};
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

    seqan3::detail::single_pass_input_view view{p};
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

    auto view = p | seqan3::views::single_pass_input | std::views::take(3);
    auto it = view.begin();

    if constexpr (std::is_same_v<std::ranges::range_value_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it, '1');
        EXPECT_EQ(*++it, '2');
        EXPECT_EQ(*++it, '3');
    }
    else
    {
        EXPECT_EQ(*it, 1);
        EXPECT_EQ(*++it, 2);
        EXPECT_EQ(*++it, 3);
    }
    ++it;
    EXPECT_TRUE(view.end() == it);
}

TYPED_TEST(single_pass_input, fn_pipeable)
{
    TypeParam p{this->data};

    auto view = p | seqan3::views::single_pass_input | std::views::take(3);
    auto it = view.begin();
    if constexpr (std::is_same_v<std::ranges::range_value_t<TypeParam>, char>)
    {
        EXPECT_EQ(*it, '1');
        EXPECT_EQ(*++it, '2');
        EXPECT_EQ(*++it, '3');
    }
    else
    {
        EXPECT_EQ(*it, 1);
        EXPECT_EQ(*++it, 2);
        EXPECT_EQ(*++it, 3);
    }
    ++it;
    EXPECT_TRUE(view.end() == it);
}
