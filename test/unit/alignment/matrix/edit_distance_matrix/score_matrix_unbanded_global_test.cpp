// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/debug_matrix.hpp>
#include <seqan3/alignment/matrix/detail/edit_distance_score_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/matrix_concept.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "../../pairwise/fixture/global_edit_distance_unbanded.hpp"
#include "edit_distance_score_matrix.hpp"

TEST(global, empty)
{
    matrix_type<false, false> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(global, epsilon)
{
    matrix_type<false, false> matrix{1u};

    matrix.add_column({}, {});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0}};

    EXPECT_EQ(result, expect);
}

TEST(global, epsilon_row)
{
    matrix_type<false, false> matrix{1u};

    matrix.add_column({}, {});
    matrix.add_column({}, {});
    matrix.add_column({}, {});
    matrix.add_column({}, {});
    matrix.add_column({}, {});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -1, -2, -3, -4}};

    EXPECT_EQ(result, expect);
}

TEST(global, single_word)
{
    matrix_type<false, false> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u}, {0b0000'0000u});
    matrix.add_column({0b1111'1110u}, {0b0000'0001u});
    matrix.add_column({0b1110'1100u}, {0b0000'0001u});
    matrix.add_column({0b1101'1100u}, {0b0010'0011u});
    matrix.add_column({0b1001'1000u}, {0b0000'0011u});
    matrix.add_column({0b1011'1000u}, {0b0100'0111u});
    matrix.add_column({0b0011'0000u}, {0b0000'0111u});
    matrix.add_column({0b0111'0000u}, {0b1000'1111u});
    matrix.add_column({0b0110'0000u}, {0b0000'1111u});
    matrix.add_column({0b1110'0000u}, {0b0001'1111u});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
                                         {-1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
                                         {-2, -1, -1, -1, -2, -3, -4, -5, -6, -7},
                                         {-3, -2, -2, -2, -2, -2, -3, -4, -5, -6},
                                         {-4, -3, -3, -3, -3, -3, -3, -3, -4, -5},
                                         {-5, -4, -3, -4, -4, -4, -4, -4, -4, -4},
                                         {-6, -5, -4, -3, -4, -5, -5, -5, -5, -5},
                                         {-7, -6, -5, -4, -4, -4, -5, -6, -6, -6},
                                         {-8, -7, -6, -5, -5, -5, -5, -5, -6, -7}};

    EXPECT_EQ(result, expect);
}

TEST(global, multiple_words)
{
    matrix_type<false, false> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1111'1110u, 0b1111'1111u, 0b1u}, {0b0000'0001u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1111'1000u, 0b1111'1111u, 0b1u}, {0b0000'0001u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1110'0000u, 0b1111'1111u, 0b1u}, {0b0000'0001u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1000'0000u, 0b1111'1111u, 0b1u}, {0b0000'0001u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0000'0100u, 0b1111'1110u, 0b1u}, {0b0000'0011u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0001'1000u, 0b1111'1000u, 0b1u}, {0b0000'0111u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0110'0000u, 0b1110'0000u, 0b1u}, {0b0000'0111u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1000'0000u, 0b1000'0001u, 0b1u}, {0b0000'0111u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0000'0000u, 0b0000'0110u, 0b0u}, {0b0000'0111u, 0b0000'0000u, 0b0u});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
                                         {-1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
                                         {-2, -1, -1, -2, -3, -3, -4, -5, -6, -7},
                                         {-3, -2, -1, -2, -3, -4, -3, -4, -5, -6},
                                         {-4, -3, -2, -2, -3, -4, -4, -4, -5, -6},
                                         {-5, -4, -3, -2, -3, -4, -5, -4, -5, -6},
                                         {-6, -5, -4, -3, -3, -4, -5, -5, -5, -6},
                                         {-7, -6, -5, -4, -3, -4, -5, -6, -5, -6},
                                         {-8, -7, -6, -5, -4, -4, -5, -6, -6, -6},
                                         {-9, -8, -7, -6, -5, -4, -5, -6, -7, -6},
                                         {-10, -9, -8, -7, -6, -5, -5, -6, -7, -7},
                                         {-11, -10, -9, -8, -7, -6, -5, -6, -7, -8},
                                         {-12, -11, -10, -9, -8, -7, -6, -6, -7, -8},
                                         {-13, -12, -11, -10, -9, -8, -7, -6, -7, -8},
                                         {-14, -13, -12, -11, -10, -9, -8, -7, -7, -8},
                                         {-15, -14, -13, -12, -11, -10, -9, -8, -7, -8},
                                         {-16, -15, -14, -13, -12, -11, -10, -9, -8, -8},
                                         {-17, -16, -15, -14, -13, -12, -11, -10, -9, -8}};

    EXPECT_EQ(result, expect);
}
