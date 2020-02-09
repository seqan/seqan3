// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/alignment/matrix/edit_distance_score_matrix_full.hpp>
#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "edit_distance_score_matrix.hpp"

TEST(max_rows, score_mismatch)
{
    using matrix_t = matrix_type<true, false>;

    // If the score mismatches max_errors, the row index obtained by score_mask and last_block contains INF.
    // -0, -1, -2,*-3,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0000'0000u, 0u, 3u, 2u), 0u);
    // -0, -1, -2, -3,
    // -1, -2, -3,*-4,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0000'0001u, 0u, 4u, 3u), 1u);
    // -1, -2, -3, -4,
    // -2, -3,*-4, -5,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0000'0010u, 0u, 4u, 3u), 2u);
    // -2, -3, -4, -5,
    // -3,*-4, -5, -6,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0000'0100u, 0u, 4u, 3u), 3u);
    // -3, -4, -5, -6,
    //*-4, -5, -6, -7,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0000'1000u, 0u, 4u, 3u), 4u);
    // -4, -5, -6, -7,
    // -5, -6, -7,*-8,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0001'0000u, 0u, 8u, 7u), 5u);
    // -5, -6, -7, -8,
    // -6, -7,*-8, -9,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0010'0000u, 0u, 8u, 7u), 6u);
    // -6, -7, -8, -9,
    // -7,*-8, -9,-10,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b0100'0000u, 0u, 8u, 7u), 7u);
    // -7, -8, -9,-10,
    //*-8, -9,-10,-11,
    // ...
    EXPECT_EQ(matrix_t::max_rows(0b1000'0000u, 0u, 8u, 7u), 8u);
    // ...
    // -8, -9,-10,-11,
    //*-9,-10,-11,-12,
    EXPECT_EQ(matrix_t::max_rows(0b0000'0001u, 1u, 9u, 8u), 9u);
}

TEST(max_rows, score_match)
{
    using matrix_t = matrix_type<true, false>;

    // If the score mismatches max_errors, the row index obtained by score_mask and last_block contains INF.
    // -0, -1, -2,*-3,
    EXPECT_EQ(matrix_t::max_rows(0b0000'0000u, 0u, 3u, 3u), 1u);
    // -0, -1, -2, -3,
    // -1,*-2, -3, -4,
    EXPECT_EQ(matrix_t::max_rows(0b0000'0001u, 0u, 2u, 4u), 2u);
    // ...
    // -1, -2, -3, -4,
    // -2, -3,*-4, -5,
    EXPECT_EQ(matrix_t::max_rows(0b0000'0010u, 0u, 4u, 4u), 3u);
    // ...
    // -2, -3, -4, -5,
    // -3,*-4, -5, -6,
    EXPECT_EQ(matrix_t::max_rows(0b0000'0100u, 0u, 4u, 4u), 4u);
    // ...
    // -3, -4, -5, -6,
    //*-4, -5, -6, -7,
    EXPECT_EQ(matrix_t::max_rows(0b0000'1000u, 0u, 4u, 4u), 5u);
    // ...
    // -4, -5, -6, -7,
    // -5, -6, -7,*-8,
    EXPECT_EQ(matrix_t::max_rows(0b0001'0000u, 0u, 8u, 10u), 6u);
    // ...
    // -5, -6, -7, -8,
    // -6, -7,*-8, -9,
    EXPECT_EQ(matrix_t::max_rows(0b0010'0000u, 0u, 8u, 10u), 7u);
    // ...
    // -6, -7, -8, -9,
    // -7,*-8, -9,-10,
    EXPECT_EQ(matrix_t::max_rows(0b0100'0000u, 0u, 8u, 10u), 8u);
    // ...
    // -7, -8, -9,-10,
    //*-8, -9,-10,-11,
    EXPECT_EQ(matrix_t::max_rows(0b1000'0000u, 0u, 8u, 8u), 9u);
    // ...
    // -8, -9,-10,-11,
    //*-9,-10,-11,-12,
    EXPECT_EQ(matrix_t::max_rows(0b0000'0001u, 1u, 9u, 9u), 10u);
}
