// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>

//------------------------------------------------------------------------------
// score cell proxy
//------------------------------------------------------------------------------

struct affine_cell_proxy_test : public ::testing::Test
{
    using score_cell_type = std::tuple<int, int const &, int &>;
    using cell_type = seqan3::detail::affine_cell_proxy<score_cell_type>;

    int best_score{4};
    int horizontal_score{-1};
    int vertical_score{10};

    cell_type affine_cell{best_score, horizontal_score, vertical_score};
};

TEST_F(affine_cell_proxy_test, construction)
{
    int lvalue_variable{8};
    cell_type other_cell{1, lvalue_variable, lvalue_variable};

    EXPECT_EQ(std::get<0>(other_cell), 1);
    EXPECT_EQ(std::get<1>(other_cell), 8);
    EXPECT_EQ(std::get<2>(other_cell), 8);

    cell_type other_cell2{affine_cell};
    EXPECT_EQ(std::get<0>(other_cell2), best_score);
    EXPECT_EQ(std::get<1>(other_cell2), horizontal_score);
    EXPECT_EQ(std::get<2>(other_cell2), vertical_score);
}

TEST_F(affine_cell_proxy_test, assignment)
{
    int lvalue_variable{8};
    seqan3::detail::affine_cell_proxy<std::tuple<int, int, int>> other_cell{1, lvalue_variable, lvalue_variable};

    EXPECT_EQ(std::get<0>(other_cell), 1);
    EXPECT_EQ(std::get<1>(other_cell), 8);
    EXPECT_EQ(std::get<2>(other_cell), 8);

    other_cell = affine_cell;
    EXPECT_EQ(std::get<0>(other_cell), best_score);
    EXPECT_EQ(std::get<1>(other_cell), horizontal_score);
    EXPECT_EQ(std::get<2>(other_cell), vertical_score);
}

TEST_F(affine_cell_proxy_test, best_score)
{
    EXPECT_EQ(affine_cell.best_score(), best_score);
    EXPECT_EQ(std::as_const(affine_cell).best_score(), best_score);
    EXPECT_EQ(std::move(affine_cell).best_score(), best_score);
    EXPECT_EQ(std::move(std::as_const(affine_cell)).best_score(), best_score);
    EXPECT_TRUE((std::same_as<decltype(affine_cell.best_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(affine_cell).best_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(affine_cell).best_score()), int &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(affine_cell)).best_score()), int const &&>));
}

TEST_F(affine_cell_proxy_test, horizontal_score)
{
    EXPECT_EQ(affine_cell.horizontal_score(), horizontal_score);
    EXPECT_EQ(std::as_const(affine_cell).horizontal_score(), horizontal_score);
    EXPECT_EQ(std::move(affine_cell).horizontal_score(), horizontal_score);
    EXPECT_EQ(std::move(std::as_const(affine_cell)).horizontal_score(), horizontal_score);
    EXPECT_TRUE((std::same_as<decltype(affine_cell.horizontal_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(affine_cell).horizontal_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(affine_cell).horizontal_score()), int const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(affine_cell)).horizontal_score()), int const &>));
}

TEST_F(affine_cell_proxy_test, vertical_score)
{
    EXPECT_EQ(affine_cell.vertical_score(), vertical_score);
    EXPECT_EQ(std::as_const(affine_cell).vertical_score(), vertical_score);
    EXPECT_EQ(std::move(affine_cell).vertical_score(), vertical_score);
    EXPECT_EQ(std::move(std::as_const(affine_cell)).vertical_score(), vertical_score);
    EXPECT_TRUE((std::same_as<decltype(affine_cell.vertical_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(affine_cell).vertical_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(affine_cell).vertical_score()), int &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(affine_cell)).vertical_score()), int &>));
}

TEST_F(affine_cell_proxy_test, tuple_size)
{
    EXPECT_EQ(std::tuple_size_v<cell_type>, 3u);
}

TEST_F(affine_cell_proxy_test, tuple_element)
{
    EXPECT_TRUE((std::same_as<std::tuple_element_t<0, cell_type>, int>));
    EXPECT_TRUE((std::same_as<std::tuple_element_t<1, cell_type>, int const &>));
    EXPECT_TRUE((std::same_as<std::tuple_element_t<2, cell_type>, int &>));
}

TEST_F(affine_cell_proxy_test, tuple_like_concept)
{
    EXPECT_TRUE(seqan3::tuple_like<cell_t>);
    EXPECT_TRUE(seqan3::detail::affine_score_cell<cell_t>);
}
