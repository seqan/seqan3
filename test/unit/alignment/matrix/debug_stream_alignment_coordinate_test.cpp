// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>

using namespace seqan3;

TEST(debug_stream_test, advanceable_alignment_coordinate)
{
    using not_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::row>;
    using col_incrementable =
        detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::column>;

    not_incrementable co_not{detail::column_index_type{10u}, detail::row_index_type{5u}};
    col_incrementable co_col{detail::column_index_type{10u}, detail::row_index_type{5u}};
    row_incrementable co_row{detail::column_index_type{10u}, detail::row_index_type{5u}};

    EXPECT_TRUE((detail::is_value_specialisation_of_v<not_incrementable, detail::advanceable_alignment_coordinate>));
    EXPECT_TRUE((detail::is_value_specialisation_of_v<col_incrementable, detail::advanceable_alignment_coordinate>));
    EXPECT_TRUE((detail::is_value_specialisation_of_v<row_incrementable, detail::advanceable_alignment_coordinate>));

    std::stringstream sstream{};
    debug_stream_type dstream{sstream};
    dstream << co_not;
    dstream << co_col;
    dstream << co_row;
    EXPECT_EQ(sstream.str(), "(10,5)(10,5)(10,5)");

    EXPECT_EQ(co_not, co_not);
    EXPECT_EQ(co_col, co_col);
    EXPECT_EQ(co_row, co_row);
}

TEST(debug_stream_test, alignment_coordinate)
{
    alignment_coordinate co_align{detail::column_index_type{10u}, detail::row_index_type{5u}};

    EXPECT_FALSE((detail::is_value_specialisation_of_v<decltype(co_align), detail::advanceable_alignment_coordinate>));

    std::stringstream sstream{};
    debug_stream_type dstream{sstream};
    dstream << co_align;
    EXPECT_EQ(sstream.str(), "(10,5)");

    EXPECT_EQ(co_align, co_align);
}
