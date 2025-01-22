// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>

TEST(debug_stream_test, advanceable_alignment_coordinate)
{
    using not_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;
    using col_incrementable = seqan3::detail::advanceable_alignment_coordinate<
        seqan3::detail::advanceable_alignment_coordinate_state::column>;

    not_incrementable co_not{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    col_incrementable co_col{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    row_incrementable co_row{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};

    EXPECT_TRUE((seqan3::detail::is_value_specialisation_of_v<not_incrementable,
                                                              seqan3::detail::advanceable_alignment_coordinate>));
    EXPECT_TRUE((seqan3::detail::is_value_specialisation_of_v<col_incrementable,
                                                              seqan3::detail::advanceable_alignment_coordinate>));
    EXPECT_TRUE((seqan3::detail::is_value_specialisation_of_v<row_incrementable,
                                                              seqan3::detail::advanceable_alignment_coordinate>));

    std::stringstream sstream{};
    seqan3::debug_stream_type dstream{sstream};
    dstream << co_not;
    dstream << co_col;
    dstream << co_row;
    EXPECT_EQ(sstream.str(), "(10,5)(10,5)(10,5)");

    EXPECT_EQ(co_not, co_not);
    EXPECT_EQ(co_col, co_col);
    EXPECT_EQ(co_row, co_row);
}
