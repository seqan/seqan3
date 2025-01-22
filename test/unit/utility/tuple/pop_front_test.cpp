// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <concepts>
#include <tuple>
#include <type_traits>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/utility/tuple/pod_tuple.hpp>
#include <seqan3/utility/tuple/pop_front.hpp>
#include <seqan3/utility/tuple/split.hpp>

#include "my_tuple.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"

struct bar : public seqan3::detail::strong_type<unsigned, bar>
{
    using seqan3::detail::strong_type<unsigned, bar>::strong_type;
};

template <typename T>
class tuple_pop_front : public ::testing::Test
{
public:
    T value;
};

using tuple_utility_types =
    ::testing::Types<std::tuple<int, long, bar, float>, seqan3::pod_tuple<int, long, bar, float>>;

TYPED_TEST_SUITE(tuple_pop_front, tuple_utility_types, );

TYPED_TEST(tuple_pop_front, lvalue)
{
    TypeParam t{1, 10l, bar{2}, 2.1};
    auto res = seqan3::tuple_pop_front(t);

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);

    auto res2 = seqan3::tuple_pop_front(seqan3::tuple_pop_front(seqan3::tuple_pop_front(res)));

    EXPECT_EQ(std::tuple_size_v<decltype(res2)>, 0u);
}

TYPED_TEST(tuple_pop_front, const_lvalue)
{
    TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
    auto res = seqan3::tuple_pop_front(t);

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);
}

TYPED_TEST(tuple_pop_front, rvalue)
{
    TypeParam t{1, 10l, bar{2}, 2.1};
    auto res = seqan3::tuple_pop_front(std::move(t));

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);
}

TYPED_TEST(tuple_pop_front, const_rvalue)
{
    TypeParam const t{TypeParam{1, 10l, bar{2}, 2.1}};
    auto res = seqan3::tuple_pop_front(std::move(t));

    EXPECT_EQ(std::tuple_size_v<decltype(res)>, 3u);

    EXPECT_EQ(std::get<0>(res), 10l);
    EXPECT_EQ(std::get<1>(res).get(), 2u);
    EXPECT_FLOAT_EQ(std::get<2>(res), 2.1);
}

TYPED_TEST(tuple_pop_front, tuple_split_and_pop)
{
    std::tuple t{float{2.1}};
    {
        auto [left, right] = seqan3::tuple_split<float>(t);

        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(left)>>, 0u);
        EXPECT_EQ(std::tuple_size_v<std::remove_reference_t<decltype(right)>>, 1u);

        using left_tuple_t = std::remove_reference_t<decltype(left)>;
        using right_tuple_t = std::remove_reference_t<decltype(seqan3::tuple_pop_front(right))>;

        using left_t = seqan3::detail::transfer_template_args_onto_t<left_tuple_t, seqan3::type_list>;
        using right_t = seqan3::detail::transfer_template_args_onto_t<right_tuple_t, seqan3::type_list>;

        EXPECT_TRUE((std::is_same_v<left_t, seqan3::type_list<>>));
        EXPECT_TRUE((std::is_same_v<right_t, seqan3::type_list<>>));

        auto v = std::tuple_cat(left, std::tuple{1}, seqan3::tuple_pop_front(right));

        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(v)>, std::tuple<int>>));
    }
}

#pragma GCC diagnostic pop
