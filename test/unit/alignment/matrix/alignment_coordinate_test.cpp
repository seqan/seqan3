// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/test/pretty_printing.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
TEST(alignment_coordinate, basic)
{
    EXPECT_TRUE(std::is_default_constructible<seqan3::alignment_coordinate>::value);
    EXPECT_TRUE(std::is_copy_constructible<seqan3::alignment_coordinate>::value);
    EXPECT_TRUE(std::is_copy_assignable<seqan3::alignment_coordinate>::value);
    EXPECT_TRUE(std::is_move_constructible<seqan3::alignment_coordinate>::value);
    EXPECT_TRUE(std::is_move_assignable<seqan3::alignment_coordinate>::value);
    EXPECT_TRUE(std::is_destructible<seqan3::alignment_coordinate>::value);

    using not_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::none>;
    using row_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::row>;
    using col_incrementable =
        seqan3::detail::advanceable_alignment_coordinate<seqan3::detail::advanceable_alignment_coordinate_state::column>;

    not_incrementable co_not{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    col_incrementable co_col{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    row_incrementable co_row{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};

    seqan3::alignment_coordinate test1{co_not};
    EXPECT_EQ(test1.first, 10u);
    EXPECT_EQ(test1.second, 5u);

    seqan3::alignment_coordinate test2{co_col};
    EXPECT_EQ(test2.first, 10u);
    EXPECT_EQ(test2.second, 5u);

    seqan3::alignment_coordinate test3{co_row};
    EXPECT_EQ(test3.first, 10u);
    EXPECT_EQ(test3.second, 5u);

    seqan3::alignment_coordinate test4{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    EXPECT_EQ(test4.first, 10u);
    EXPECT_EQ(test4.second, 5u);
}

TEST(alignment_coordinate, matrix_coordainte_conversion)
{
    seqan3::alignment_coordinate co{seqan3::detail::column_index_type{10u}, seqan3::detail::row_index_type{5u}};
    seqan3::detail::matrix_coordinate mc = co;

    EXPECT_EQ(mc.col, 10u);
    EXPECT_EQ(mc.row, 5u);
}
#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
