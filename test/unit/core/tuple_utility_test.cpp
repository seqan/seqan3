// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <tuple>

#include "my_tuple.hpp"

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/std/concepts>

using namespace seqan3;

struct bar : public detail::strong_type<unsigned, bar>
{
    using detail::strong_type<unsigned, bar>::strong_type;
};

template <typename T>
class tuple_utility : public ::testing::Test
{
public:

    T value;
};

using tuple_utility_types = ::testing::Types<std::tuple<int, long, bar, float>,
                                             pod_tuple<int, long, bar, float>>;

TYPED_TEST_CASE(tuple_utility, tuple_utility_types);

TYPED_TEST(tuple_utility, tuple_type_list)
{
    {
        using list = typename detail::tuple_type_list<my_tuple>::type;
        EXPECT_TRUE((std::is_same_v<list, type_list<int, float>>));
    }

    {
        using list = detail::tuple_type_list_t<TypeParam>;
        EXPECT_TRUE((std::is_same_v<list, type_list<int, long, bar, float>>));
    }
}

TYPED_TEST(tuple_utility, TupleLike)
{
    EXPECT_TRUE(TupleLike<TypeParam>);
    EXPECT_TRUE(TupleLike<std::tuple<>>);
    EXPECT_TRUE(TupleLike<my_tuple>);
    EXPECT_FALSE(TupleLike<int>);
}

TYPED_TEST(tuple_utility, detail_split)
{
    TypeParam t{1, 10l, bar{2}, 2.1};

    {
        auto res = detail::tuple_split<0>(t, std::make_index_sequence<0>{});
        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 0u);
    }

    {
        auto res = detail::tuple_split<2>(t, std::make_index_sequence<2>{});
        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, decltype(res)>, bar>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, decltype(res)>, float>));
        EXPECT_EQ(std::get<0>(res).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<1>(res), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_pos_lvalue)
{
    TypeParam t{1, 10l, bar{2}, 2.1};
    {
        auto res = tuple_split<0>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }

    {
        auto res = tuple_split<1>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 1u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 3u);

        EXPECT_EQ(std::get<0>(std::get<0>(res)), 1);
        EXPECT_EQ(std::get<0>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<1>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<2>(std::get<1>(res)), 2.1);
    }

    {
        auto res = tuple_split<3>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 3u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 1u);
    }

    {
        auto res = tuple_split<4>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 4u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 0u);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_pos_const_lvalue)
{
    TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
    EXPECT_TRUE((std::is_same_v<decltype(t), TypeParam const>));
    {
        auto res = tuple_split<0>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_pos_rvalue)
{
    {
        auto res = tuple_split<0>(TypeParam{1, 10l, bar{2}, 2.1});

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_pos_const_rvalue)
{
    {
        TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
        EXPECT_TRUE((std::is_same_v<decltype(t), TypeParam const>));
        auto res = tuple_split<0>(std::move(t));

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_type_lvalue)
{
    TypeParam t{1, 10l, bar{2}, 2.1};
    {
        auto res = tuple_split<int>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }

    {
        auto res = tuple_split<long int>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 1u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 3u);

        EXPECT_EQ(std::get<0>(std::get<0>(res)), 1);
        EXPECT_EQ(std::get<0>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<1>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<2>(std::get<1>(res)), 2.1);
    }

    {
        auto res = tuple_split<float>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 3u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 1u);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_type_const_lvalue)
{
    TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
    EXPECT_TRUE((std::is_same_v<decltype(t), TypeParam const>));
    {
        auto res = tuple_split<int>(t);

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_type_rvalue)
{
    {
        auto res = tuple_split<int>(TypeParam{1, 10l, bar{2}, 2.1});

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_split_by_type_const_rvalue)
{
    {
        TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
        EXPECT_TRUE((std::is_same_v<decltype(t), TypeParam const>));
        auto res = tuple_split<int>(std::move(t));

        EXPECT_EQ(std::tuple_size_v<decltype(res)>, 2u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<0>(res))>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(std::get<1>(res))>>, 4u);

        EXPECT_EQ(std::get<0>(std::get<1>(res)), 1);
        EXPECT_EQ(std::get<1>(std::get<1>(res)), 10l);
        EXPECT_EQ(std::get<2>(std::get<1>(res)).get(), 2u);
        EXPECT_FLOAT_EQ(std::get<3>(std::get<1>(res)), 2.1);
    }
}

TYPED_TEST(tuple_utility, tuple_pop_front_lvalue)
{
    TypeParam t{1, 10l, bar{2}, 2.1};
    auto res = tuple_pop_front(t);

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);

    auto res2 = tuple_pop_front(tuple_pop_front(tuple_pop_front(res)));

    EXPECT_EQ(std::tuple_size_v<decltype(res2)>, 0u);
}

TYPED_TEST(tuple_utility, tuple_pop_front_const_lvalue)
{
    TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
    auto res = tuple_pop_front(t);

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);
}

TYPED_TEST(tuple_utility, tuple_pop_front_rvalue)
{
    TypeParam t{1, 10l, bar{2}, 2.1};
    auto res = tuple_pop_front(std::move(t));

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);
}

TYPED_TEST(tuple_utility, tuple_pop_front_const_rvalue)
{
    TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
    auto res = tuple_pop_front(std::move(t));

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);
}

TYPED_TEST(tuple_utility, tuple_split_and_pop)
{
    std::tuple t{float{2.1}};
    {
        auto [left, right] = tuple_split<float>(t);

        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(left)>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(right)>>, 1u);

        using left_t = detail::transfer_template_args_onto_t<remove_cvref_t<decltype(left)>, type_list>;
        using right_t = detail::transfer_template_args_onto_t<remove_cvref_t<decltype(tuple_pop_front(right))>, type_list>;

        EXPECT_TRUE((std::is_same_v<left_t, type_list<>>));
        EXPECT_TRUE((std::is_same_v<right_t, type_list<>>));

        auto v = std::tuple_cat(left, std::tuple{1}, tuple_pop_front(right));

        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(v)>, std::tuple<int>>));
    }
}
