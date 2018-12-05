// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

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

TYPED_TEST(tuple_utility, tuple_like_concept)
{
    EXPECT_TRUE(tuple_like_concept<TypeParam>);
    EXPECT_TRUE(tuple_like_concept<std::tuple<>>);
    EXPECT_TRUE(tuple_like_concept<my_tuple>);
    EXPECT_FALSE(tuple_like_concept<int>);
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
