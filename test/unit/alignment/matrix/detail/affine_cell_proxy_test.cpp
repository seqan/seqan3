// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>

struct affine_cell_proxy_test : public ::testing::Test
{
    using base_t = std::tuple<int, int const &, int &>;
    using cell_t = seqan3::detail::affine_cell_proxy<base_t>;

    int optimal{4};
    int horizontal{-1};
    int vertical{10};

    cell_t affine_cell{optimal, horizontal, vertical};
};

TEST_F(affine_cell_proxy_test, construction)
{
    int lvalue_variable{8};
    cell_t other_cell{1, lvalue_variable, lvalue_variable};

    EXPECT_EQ(std::get<0>(other_cell), 1);
    EXPECT_EQ(std::get<1>(other_cell), 8);
    EXPECT_EQ(std::get<2>(other_cell), 8);

    cell_t other_cell2{affine_cell};
    EXPECT_EQ(std::get<0>(other_cell2), 4);
    EXPECT_EQ(std::get<1>(other_cell2), -1);
    EXPECT_EQ(std::get<2>(other_cell2), 10);
}

TEST_F(affine_cell_proxy_test, assignment)
{
    int lvalue_variable{8};
    seqan3::detail::affine_cell_proxy<std::tuple<int, int, int>> other_cell{1, lvalue_variable, lvalue_variable};

    EXPECT_EQ(std::get<0>(other_cell), 1);
    EXPECT_EQ(std::get<1>(other_cell), 8);
    EXPECT_EQ(std::get<2>(other_cell), 8);

    other_cell = affine_cell;
    EXPECT_EQ(std::get<0>(other_cell), 4);
    EXPECT_EQ(std::get<1>(other_cell), -1);
    EXPECT_EQ(std::get<2>(other_cell), 10);
}

TEST_F(affine_cell_proxy_test, best_score)
{
    EXPECT_EQ(affine_cell.best_score(), 4);
    EXPECT_EQ(std::as_const(affine_cell).best_score(), 4);
    EXPECT_EQ(std::move(affine_cell).best_score(), 4);
    EXPECT_EQ(std::move(std::as_const(affine_cell)).best_score(), 4);
    EXPECT_TRUE((std::same_as<decltype(affine_cell.best_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(affine_cell).best_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(affine_cell).best_score()), int &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(affine_cell)).best_score()), int const &&>));
}

TEST_F(affine_cell_proxy_test, horizontal_score)
{
    EXPECT_EQ(affine_cell.horizontal_score(), -1);
    EXPECT_EQ(std::as_const(affine_cell).horizontal_score(), -1);
    EXPECT_EQ(std::move(affine_cell).horizontal_score(), -1);
    EXPECT_EQ(std::move(std::as_const(affine_cell)).horizontal_score(), -1);
    EXPECT_TRUE((std::same_as<decltype(affine_cell.horizontal_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(affine_cell).horizontal_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(affine_cell).horizontal_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(affine_cell)).horizontal_score()), int const &>));
}

TEST_F(affine_cell_proxy_test, vertical_score)
{
    EXPECT_EQ(affine_cell.vertical_score(), 10);
    EXPECT_EQ(std::as_const(affine_cell).vertical_score(), 10);
    EXPECT_EQ(std::move(affine_cell).vertical_score(), 10);
    EXPECT_EQ(std::move(std::as_const(affine_cell)).vertical_score(), 10);
    EXPECT_TRUE((std::same_as<decltype(affine_cell.vertical_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(affine_cell).vertical_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(affine_cell).vertical_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(affine_cell)).vertical_score()), int &>));
}

TEST_F(affine_cell_proxy_test, tuple_size)
{
    EXPECT_EQ(std::tuple_size_v<cell_t>, 3u);
}

TEST_F(affine_cell_proxy_test, tuple_element)
{
    EXPECT_TRUE((std::same_as<std::tuple_element_t<0, cell_t>, int>));
    EXPECT_TRUE((std::same_as<std::tuple_element_t<1, cell_t>, int const &>));
    EXPECT_TRUE((std::same_as<std::tuple_element_t<2, cell_t>, int &>));
}

TEST_F(affine_cell_proxy_test, tuple_like_concept)
{
    EXPECT_TRUE(seqan3::tuple_like<cell_t>);
}
