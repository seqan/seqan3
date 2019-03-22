// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/alignment_optimum.hpp>

using namespace seqan3;

TEST(alignment_optimum, construction)
{
    EXPECT_TRUE((std::is_default_constructible<detail::alignment_optimum<int32_t>>::value));
    EXPECT_TRUE((std::is_copy_constructible<detail::alignment_optimum<int32_t>>::value));
    EXPECT_TRUE((std::is_copy_assignable<detail::alignment_optimum<int32_t>>::value));
    EXPECT_TRUE((std::is_move_constructible<detail::alignment_optimum<int32_t>>::value));
    EXPECT_TRUE((std::is_move_assignable<detail::alignment_optimum<int32_t>>::value));
    EXPECT_TRUE((std::is_destructible<detail::alignment_optimum<int32_t>>::value));
}

TEST(alignment_optimum, type_deduction)
{
    detail::alignment_optimum def_ao{};
    EXPECT_TRUE((std::is_same_v<decltype(def_ao), detail::alignment_optimum<int32_t>>));

    alignment_coordinate coordinate{};
    coordinate.first = 2u;
    coordinate.second = 3u;
    detail::alignment_optimum ao{10, coordinate};
    EXPECT_TRUE((std::is_same_v<decltype(ao), detail::alignment_optimum<int32_t>>));
}

TEST(alignment_optimum, access)
{
    detail::alignment_optimum def_ao{};

    EXPECT_EQ(def_ao.score, std::numeric_limits<int32_t>::lowest());
    EXPECT_EQ(static_cast<size_t>(def_ao.coordinate.first), static_cast<size_t>(0));
    EXPECT_EQ(static_cast<size_t>(def_ao.coordinate.second), static_cast<size_t>(0));

    alignment_coordinate coordinate{};
    coordinate.first = 2u;
    coordinate.second = 3u;
    detail::alignment_optimum ao{10, coordinate};
    EXPECT_EQ(ao.score, 10);
    EXPECT_EQ(static_cast<size_t>(ao.coordinate.first), static_cast<size_t>(2));
    EXPECT_EQ(static_cast<size_t>(ao.coordinate.second), static_cast<size_t>(3));
}

TEST(alignment_optimum, max)
{
    detail::alignment_optimum def_ao{};

    EXPECT_EQ(def_ao.score, std::numeric_limits<int32_t>::lowest());
    EXPECT_EQ(static_cast<size_t>(def_ao.coordinate.first), static_cast<size_t>(0));
    EXPECT_EQ(static_cast<size_t>(def_ao.coordinate.second), static_cast<size_t>(0));

    alignment_coordinate coordinate{};
    coordinate.first = 2u;
    coordinate.second = 3u;
    detail::alignment_optimum ao{10, coordinate};
    EXPECT_EQ(ao.score, 10);
    EXPECT_EQ(static_cast<size_t>(ao.coordinate.first), static_cast<size_t>(2));
    EXPECT_EQ(static_cast<size_t>(ao.coordinate.second), static_cast<size_t>(3));

    auto val = std::max(def_ao, ao, detail::alignment_optimum_compare_less{});
    EXPECT_EQ(ao.score, 10);
    EXPECT_EQ(static_cast<size_t>(val.coordinate.first), static_cast<size_t>(2));
    EXPECT_EQ(static_cast<size_t>(val.coordinate.second), static_cast<size_t>(3));
}
