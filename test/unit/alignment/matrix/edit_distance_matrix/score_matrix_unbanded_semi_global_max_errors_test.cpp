// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/debug_matrix.hpp>
#include <seqan3/alignment/matrix/detail/edit_distance_score_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/matrix_concept.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "edit_distance_score_matrix.hpp"

TEST(semi_global_max_errors, empty)
{
    matrix_type<true, true> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, epsilon)
{
    matrix_type<true, true> matrix{1u};

    matrix.add_column({}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, epsilon_row)
{
    matrix_type<true, true> matrix{1u};

    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -0, -0, -0, -0}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, single_word)
{
    matrix_type<true, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u}, {0b0000'0000u}, 6u);
    matrix.add_column({0b1111'1110u}, {0b0000'0000u}, 7u);
    matrix.add_column({0b1110'1110u}, {0b0000'0000u}, 8u);
    matrix.add_column({0b1101'1101u}, {0b0000'0010u}, 9u);
    matrix.add_column({0b1101'1001u}, {0b0000'0000u}, 9u);
    matrix.add_column({0b1011'1011u}, {0b0100'0100u}, 9u);
    matrix.add_column({0b0011'0011u}, {0b0000'0000u}, 9u);
    matrix.add_column({0b0111'0111u}, {0b1000'1000u}, 9u);
    matrix.add_column({0b0110'0111u}, {0b0000'0000u}, 9u);
    matrix.add_column({0b1110'1110u}, {0b0000'0000u}, 8u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -0, -0, -0, -0, -0, -0, -0, -0, -0},
                                         {-1, -0, -0, -1, -1, -1, -1, -1, -1, -0},
                                         {-2, -1, -1, -0, -1, -2, -2, -2, -2, -1},
                                         {-3, -2, -2, -1, -1, -1, -2, -3, -3, -2},
                                         {-4, -3, -3, -2, -2, -2, -2, -2, -3, -3},
                                         {-5, -4, -3, -3, -3, -3, -3, -3, -3, -3},
                                         {INF, -5, -4, -3, -3, -4, -4, -4, -4, -4},
                                         {INF, INF, -5, -4, -4, -3, -4, -5, -5, -5},
                                         {INF, INF, INF, -5, -5, -4, -4, -4, -5, INF}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global_max_errors, multiple_words)
{
    matrix_type<true, true> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 9u);
    matrix.add_column({0b1111'1110u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 10u);
    matrix.add_column({0b1111'1001u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 11u);
    matrix.add_column({0b1110'0011u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 12u);
    matrix.add_column({0b1000'0111u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 13u);
    matrix.add_column({0b0001'1110u, 0b1111'1110u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 14u);
    matrix.add_column({0b0111'1101u, 0b1111'1000u, 0b1u}, {0b0000'0010u, 0b0000'0000u, 0b0u}, 15u);
    matrix.add_column({0b1111'0001u, 0b1110'0001u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 16u);
    matrix.add_column({0b1100'0011u, 0b1000'0111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u}, 17u);
    matrix.add_column({0b0100'1110u, 0b0001'1111u, 0b0u}, {0b0001'0000u, 0b0000'0000u, 0b0u}, 18u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -0, -0, -0, -0, -0, -0, -0, -0, -0},
                                         {-1, -0, -1, -1, -1, -0, -1, -1, -1, -0},
                                         {-2, -1, -1, -2, -2, -1, -0, -1, -2, -1},
                                         {-3, -2, -1, -2, -3, -2, -1, -1, -2, -2},
                                         {-4, -3, -2, -2, -3, -3, -2, -1, -2, -3},
                                         {-5, -4, -3, -2, -3, -4, -3, -2, -2, -2},
                                         {-6, -5, -4, -3, -3, -4, -4, -3, -2, -2},
                                         {-7, -6, -5, -4, -3, -4, -5, -4, -3, -3},
                                         {-8, -7, -6, -5, -4, -4, -5, -5, -4, -3},
                                         {INF, -8, -7, -6, -5, -4, -5, -6, -5, -4},
                                         {INF, INF, -8, -7, -6, -5, -5, -6, -6, -5},
                                         {INF, INF, INF, -8, -7, -6, -5, -6, -7, -6},
                                         {INF, INF, INF, INF, -8, -7, -6, -6, -7, -7},
                                         {INF, INF, INF, INF, INF, -8, -7, -6, -7, -8},
                                         {INF, INF, INF, INF, INF, INF, -8, -7, -7, -8},
                                         {INF, INF, INF, INF, INF, INF, INF, -8, -7, -8},
                                         {INF, INF, INF, INF, INF, INF, INF, INF, -8, -8},
                                         {INF, INF, INF, INF, INF, INF, INF, INF, INF, -8}};

    EXPECT_EQ(result, expect);
}
