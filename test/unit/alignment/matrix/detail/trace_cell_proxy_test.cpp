// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <seqan3/alignment/matrix/detail/trace_cell_proxy.hpp>

using trace_t = seqan3::detail::trace_directions;

struct trace_cell_proxy_test : public ::testing::Test
{
    using base_t = std::tuple<trace_t, trace_t const &, trace_t &>;
    using cell_t = seqan3::detail::trace_cell_proxy<base_t>;

    seqan3::detail::trace_directions none = seqan3::detail::trace_directions::none;
    seqan3::detail::trace_directions diagonal = seqan3::detail::trace_directions::diagonal;
    seqan3::detail::trace_directions up = seqan3::detail::trace_directions::up;
    seqan3::detail::trace_directions left = seqan3::detail::trace_directions::left;

    cell_t trace_cell{diagonal, left, up};
};

TEST_F(trace_cell_proxy_test, construction)
{
    trace_t lvalue_variable{none};
    cell_t other_cell{seqan3::detail::trace_directions::up_open, lvalue_variable, lvalue_variable};

    EXPECT_TRUE(std::get<0>(other_cell) == seqan3::detail::trace_directions::up_open);
    EXPECT_TRUE(std::get<1>(other_cell) == none);
    EXPECT_TRUE(std::get<2>(other_cell) == none);

    cell_t other_cell2{trace_cell};
    EXPECT_TRUE(std::get<0>(other_cell2) == diagonal);
    EXPECT_TRUE(std::get<1>(other_cell2) == left);
    EXPECT_TRUE(std::get<2>(other_cell2) == up);
}

TEST_F(trace_cell_proxy_test, assignment)
{
    trace_t lvalue_variable{none};
    using tuple_t = std::tuple<trace_t, trace_t, trace_t>;
    seqan3::detail::trace_cell_proxy<tuple_t> other_cell{seqan3::detail::trace_directions::up_open,
                                                         lvalue_variable,
                                                         lvalue_variable};

    EXPECT_TRUE(std::get<0>(other_cell) == seqan3::detail::trace_directions::up_open);
    EXPECT_TRUE(std::get<1>(other_cell) == none);
    EXPECT_TRUE(std::get<2>(other_cell) == none);

    other_cell = trace_cell;
    EXPECT_TRUE(std::get<0>(other_cell) == diagonal);
    EXPECT_TRUE(std::get<1>(other_cell) == left);
    EXPECT_TRUE(std::get<2>(other_cell) == up);
}

TEST_F(trace_cell_proxy_test, trace)
{
    EXPECT_TRUE(trace_cell.trace() == diagonal);
    EXPECT_TRUE(std::as_const(trace_cell).trace() == diagonal);
    EXPECT_TRUE(cell_t{trace_cell}.trace() == diagonal);
    EXPECT_TRUE(std::move(std::as_const(trace_cell)).trace() == diagonal);
    EXPECT_TRUE((std::same_as<decltype(trace_cell.trace()), trace_t &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(trace_cell).trace()), trace_t const &>));
    EXPECT_TRUE((std::same_as<decltype(cell_t{trace_cell}.trace()), trace_t &&>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(trace_cell)).trace()), trace_t const &&>));
}

TEST_F(trace_cell_proxy_test, horizontal_trace)
{
    EXPECT_TRUE(trace_cell.horizontal_trace() == left);
    EXPECT_TRUE(std::as_const(trace_cell).horizontal_trace() == left);
    EXPECT_TRUE(cell_t{trace_cell}.horizontal_trace() == left);
    EXPECT_TRUE(std::move(std::as_const(trace_cell)).horizontal_trace() == left);
    EXPECT_TRUE((std::same_as<decltype(trace_cell.horizontal_trace()), trace_t const &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(trace_cell).horizontal_trace()), trace_t const &>));
    EXPECT_TRUE((std::same_as<decltype(cell_t{trace_cell}.horizontal_trace()), trace_t const &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(trace_cell)).horizontal_trace()), trace_t const &>));
}

TEST_F(trace_cell_proxy_test, vertical_trace)
{
    EXPECT_TRUE(trace_cell.vertical_trace() == up);
    EXPECT_TRUE(std::as_const(trace_cell).vertical_trace() == up);
    EXPECT_TRUE(cell_t{trace_cell}.vertical_trace() == up);
    EXPECT_TRUE(std::move(std::as_const(trace_cell)).vertical_trace() == up);
    EXPECT_TRUE((std::same_as<decltype(trace_cell.vertical_trace()), trace_t &>));
    EXPECT_TRUE((std::same_as<decltype(std::as_const(trace_cell).vertical_trace()), trace_t &>));
    EXPECT_TRUE((std::same_as<decltype(cell_t{trace_cell}.vertical_trace()), trace_t &>));
    EXPECT_TRUE((std::same_as<decltype(std::move(std::as_const(trace_cell)).vertical_trace()), trace_t &>));
}

TEST_F(trace_cell_proxy_test, tuple_size)
{
    EXPECT_EQ(std::tuple_size_v<cell_t>, 3u);
}

TEST_F(trace_cell_proxy_test, tuple_element)
{
    EXPECT_TRUE((std::same_as<std::tuple_element_t<0, cell_t>, trace_t>));
    EXPECT_TRUE((std::same_as<std::tuple_element_t<1, cell_t>, trace_t const &>));
    EXPECT_TRUE((std::same_as<std::tuple_element_t<2, cell_t>, trace_t &>));
}

TEST_F(trace_cell_proxy_test, tuple_like_concept)
{
    EXPECT_TRUE(seqan3::tuple_like<cell_t>);
}
