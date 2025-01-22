// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

struct A
{
    using type = int;
};

struct B;

struct C
{};

struct D
{
    static constexpr int type = 6;
};

TEST(transformation_trait_or, transformation_trait_or)
{
    using a_type = seqan3::detail::transformation_trait_or_t<A, void>;
    using b_transformation_trait_or = seqan3::detail::transformation_trait_or_t<B, void>;
    using c_transformation_trait_or = seqan3::detail::transformation_trait_or_t<C, double>;
    using d_transformation_trait_or = seqan3::detail::transformation_trait_or<D, B>::type;

    EXPECT_TRUE((std::is_same_v<a_type, int>));
    EXPECT_TRUE((std::is_same_v<b_transformation_trait_or, void>));
    EXPECT_TRUE((std::is_same_v<c_transformation_trait_or, double>));
    EXPECT_TRUE((std::is_same_v<d_transformation_trait_or, B>));
}
