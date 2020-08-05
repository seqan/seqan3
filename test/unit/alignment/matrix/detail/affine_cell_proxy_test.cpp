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
    EXPECT_TRUE(seqan3::tuple_like<cell_type>);
    EXPECT_TRUE(seqan3::detail::affine_score_cell<cell_type>);
}

//------------------------------------------------------------------------------
// combined score and trace cell proxy
//------------------------------------------------------------------------------

struct trace_cell_proxy_test : public affine_cell_proxy_test
{
    using trace_type = seqan3::detail::trace_directions;
    using trace_cell_type = std::tuple<trace_type, trace_type const &, trace_type &>;
    using cell_type = seqan3::detail::affine_cell_proxy<std::pair<score_cell_type, trace_cell_type>>;

    trace_type best_trace{trace_type::diagonal};
    trace_type horizontal_trace{trace_type::left};
    trace_type vertical_trace{trace_type::up};

    cell_type trace_cell{{best_score, horizontal_score, vertical_score},
                         {best_trace, horizontal_trace, vertical_trace}};
};

TEST_F(trace_cell_proxy_test, construction)
{
    trace_type local_trace{trace_type::none};
    int local_score{0};
    cell_type other_cell{{0, 0, local_score}, {trace_type::up_open, local_trace, local_trace}};

    EXPECT_TRUE(std::get<0>(std::get<1>(other_cell)) == trace_type::up_open);
    EXPECT_TRUE(std::get<1>(std::get<1>(other_cell)) == trace_type::none);
    EXPECT_TRUE(std::get<2>(std::get<1>(other_cell)) == trace_type::none);

    cell_type other_cell2{trace_cell};
    EXPECT_TRUE(std::get<0>(std::get<1>(other_cell2)) == best_trace);
    EXPECT_TRUE(std::get<1>(std::get<1>(other_cell2)) == horizontal_trace);
    EXPECT_TRUE(std::get<2>(std::get<1>(other_cell2)) == vertical_trace);
}

TEST_F(trace_cell_proxy_test, assignment)
{
    using local_score_cell_t = std::tuple<int, int, int>;
    using local_trace_cell_t = std::tuple<trace_type, trace_type, trace_type>;
    using other_cell_proxy_t = seqan3::detail::affine_cell_proxy<std::pair<local_score_cell_t, local_trace_cell_t>>;

    trace_type local_trace{trace_type::none};
    other_cell_proxy_t other_cell{{0, 1, 2}, {trace_type::up_open, local_trace, local_trace}};

    EXPECT_TRUE(std::get<0>(std::get<0>(other_cell)) == 0);
    EXPECT_TRUE(std::get<1>(std::get<0>(other_cell)) == 1);
    EXPECT_TRUE(std::get<2>(std::get<0>(other_cell)) == 2);

    EXPECT_TRUE(std::get<0>(std::get<1>(other_cell)) == trace_type::up_open);
    EXPECT_TRUE(std::get<1>(std::get<1>(other_cell)) == trace_type::none);
    EXPECT_TRUE(std::get<2>(std::get<1>(other_cell)) == trace_type::none);

    other_cell = trace_cell;
    EXPECT_TRUE(std::get<0>(std::get<0>(other_cell)) == best_score);
    EXPECT_TRUE(std::get<1>(std::get<0>(other_cell)) == horizontal_score);
    EXPECT_TRUE(std::get<2>(std::get<0>(other_cell)) == vertical_score);

    EXPECT_TRUE(std::get<0>(std::get<1>(other_cell)) == best_trace);
    EXPECT_TRUE(std::get<1>(std::get<1>(other_cell)) == horizontal_trace);
    EXPECT_TRUE(std::get<2>(std::get<1>(other_cell)) == vertical_trace);
}

TEST_F(trace_cell_proxy_test, best_trace)
{
    EXPECT_TRUE(trace_cell.best_trace() == best_trace);
    EXPECT_TRUE(std::as_const(trace_cell).best_trace() == best_trace);
    EXPECT_TRUE(cell_type{trace_cell}.best_trace() == best_trace);
    EXPECT_TRUE(std::move(std::as_const(trace_cell)).best_trace() == best_trace);
    EXPECT_TRUE((std::same_as<decltype(trace_cell.best_trace()), trace_type &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(trace_cell).best_trace()), trace_type const &>));
    EXPECT_TRUE((std::same_as<decltype(cell_type{trace_cell}.best_trace()), trace_type &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(trace_cell)).best_trace()), trace_type const &&>));
}

TEST_F(trace_cell_proxy_test, horizontal_trace)
{
    EXPECT_TRUE(trace_cell.horizontal_trace() == horizontal_trace);
    EXPECT_TRUE(std::as_const(trace_cell).horizontal_trace() == horizontal_trace);
    EXPECT_TRUE(cell_type{trace_cell}.horizontal_trace() == horizontal_trace);
    EXPECT_TRUE(std::move(std::as_const(trace_cell)).horizontal_trace() == horizontal_trace);
    EXPECT_TRUE((std::same_as<decltype(trace_cell.horizontal_trace()), trace_type const &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(trace_cell).horizontal_trace()), trace_type const &>));
    EXPECT_TRUE((std::same_as<decltype(cell_type{trace_cell}.horizontal_trace()), trace_type const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(trace_cell)).horizontal_trace()), trace_type const &>));
}

TEST_F(trace_cell_proxy_test, vertical_trace)
{
    EXPECT_TRUE(trace_cell.vertical_trace() == vertical_trace);
    EXPECT_TRUE(std::as_const(trace_cell).vertical_trace() == vertical_trace);
    EXPECT_TRUE(cell_type{trace_cell}.vertical_trace() == vertical_trace);
    EXPECT_TRUE(std::move(std::as_const(trace_cell)).vertical_trace() == vertical_trace);
    EXPECT_TRUE((std::same_as<decltype(trace_cell.vertical_trace()), trace_type &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(trace_cell).vertical_trace()), trace_type &>));
    EXPECT_TRUE((std::same_as<decltype(cell_type{trace_cell}.vertical_trace()), trace_type &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(trace_cell)).vertical_trace()), trace_type &>));
}
