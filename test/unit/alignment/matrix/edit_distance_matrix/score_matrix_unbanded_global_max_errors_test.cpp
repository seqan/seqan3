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

TEST(global_max_errors, empty)
{
    matrix_type<false, true> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, epsilon)
{
    matrix_type<false, true> matrix{1u};

    matrix.add_column({}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, epsilon_row)
{
    matrix_type<false, true> matrix{1u};

    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 1u);
    matrix.add_column({}, {}, 0u);
    matrix.add_column({}, {}, 0u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -1, -2,INF,INF}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_1)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u}, {0b0000'0000u}, 6u);
    matrix.add_column({0b1111'1110u}, {0b0000'0001u}, 7u);
    matrix.add_column({0b1110'1100u}, {0b0000'0001u}, 8u);
    matrix.add_column({0b1101'1100u}, {0b0010'0011u}, 9u);
    matrix.add_column({0b1001'1000u}, {0b0000'0011u}, 9u);
    matrix.add_column({0b1011'1000u}, {0b0100'0111u}, 9u);
    matrix.add_column({0b0011'0000u}, {0b0000'0111u}, 9u);
    matrix.add_column({0b0111'0000u}, {0b1000'1111u}, 9u);
    matrix.add_column({0b0110'0000u}, {0b0000'1111u}, 7u);
    matrix.add_column({0b1110'0000u}, {0b0001'1111u}, 7u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        { -0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
        { -1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
        { -2, -1, -1, -1, -2, -3, -4, -5, -6, -7},
        { -3, -2, -2, -2, -2, -2, -3, -4, -5, -6},
        { -4, -3, -3, -3, -3, -3, -3, -3, -4, -5},
        { -5, -4, -3, -4, -4, -4, -4, -4, -4, -4},
        {INF, -5, -4, -3, -4, -5, -5, -5, -5, -5},
        {INF,INF, -5, -4, -4, -4, -5, -6,INF,INF},
        {INF,INF,INF, -5, -5, -5, -5, -5,INF,INF}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_2)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u}, {0b0000'0000u}, 5u);
    matrix.add_column({0b1111'1110u}, {0b0000'0001u}, 6u);
    matrix.add_column({0b1110'1100u}, {0b0000'0001u}, 7u);
    matrix.add_column({0b1101'1100u}, {0b0010'0011u}, 8u);
    matrix.add_column({0b1001'1000u}, {0b0000'0011u}, 8u);
    matrix.add_column({0b1011'1000u}, {0b0100'0111u}, 8u);
    matrix.add_column({0b0011'0000u}, {0b0000'0111u}, 6u);
    matrix.add_column({0b0111'0000u}, {0b1000'1111u}, 6u);
    matrix.add_column({0b0110'0000u}, {0b0000'1111u}, 6u);
    matrix.add_column({0b1110'0000u}, {0b0001'1111u}, 6u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        { -0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
        { -1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
        { -2, -1, -1, -1, -2, -3, -4, -5, -6, -7},
        { -3, -2, -2, -2, -2, -2, -3, -4, -5, -6},
        { -4, -3, -3, -3, -3, -3, -3, -3, -4, -5},
        {INF, -4, -3, -4, -4, -4, -4, -4, -4, -4},
        {INF,INF, -4, -3, -4, -5,INF,INF,INF,INF},
        {INF,INF,INF, -4, -4, -4,INF,INF,INF,INF},
        {INF,INF,INF,INF,INF,INF,INF,INF,INF,INF}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_3)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    // Note that score_mask = 0b0000'1000u means that only the number of bits up until the 1 (from right-to-left) is
    // relevant. That means only the first 4 bits (from right-to-left) in 0b1010'1111u are relevant and so the X's of
    // 0bXXXX'1111u can be filled with anything. Furthermore, note that we filled "random" bits in these test cases.
    matrix.add_column({0b1010'1111u}, {0b0101'0000u}, 4u);
    matrix.add_column({0b0101'1110u}, {0b1010'0001u}, 5u);
    matrix.add_column({0b1010'1100u}, {0b0100'0001u}, 6u);
    matrix.add_column({0b0101'1100u}, {0b1010'0011u}, 7u);
    matrix.add_column({0b0101'1000u}, {0b1010'0011u}, 5u);
    matrix.add_column({0b0101'1000u}, {0b1010'0111u}, 5u);
    matrix.add_column({0b0101'0000u}, {0b1010'0111u}, 5u);
    matrix.add_column({0b0101'0000u}, {0b1010'1111u}, 5u);
    matrix.add_column({0b1010'1010u}, {0b0101'0101u}, 0u);
    matrix.add_column({0b1010'1010u}, {0b0101'0101u}, 0u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        { -0, -1, -2, -3, -4, -5, -6, -7,INF,INF},
        { -1, -0, -1, -2, -3, -4, -5, -6,INF,INF},
        { -2, -1, -1, -1, -2, -3, -4, -5,INF,INF},
        { -3, -2, -2, -2, -2, -2, -3, -4,INF,INF},
        {INF, -3, -3, -3, -3, -3, -3, -3,INF,INF},
        {INF,INF, -3, -4,INF,INF,INF,INF,INF,INF},
        {INF,INF,INF, -3,INF,INF,INF,INF,INF,INF},
        {INF,INF,INF,INF,INF,INF,INF,INF,INF,INF},
        {INF,INF,INF,INF,INF,INF,INF,INF,INF,INF}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, multiple_words_1)
{
    matrix_type<false, true> matrix{10u};
    matrix.reserve(10u);

    matrix.add_column({0b0111'1111u}, {0b1000'0000u}, 6u);
    matrix.add_column({0b1111'1110u}, {0b0000'0001u}, 7u);
    matrix.add_column({0b1110'1100u}, {0b0000'0001u}, 8u);
    matrix.add_column({0b1101'1100u, 0b1u}, {0b0010'0011u, 0b0u}, 9u);
    matrix.add_column({0b1001'1000u, 0b1u}, {0b0000'0011u, 0b0u}, 9u);
    matrix.add_column({0b1011'1000u, 0b1u}, {0b0100'0111u, 0b0u}, 9u);
    matrix.add_column({0b0011'0000u, 0b1u}, {0b0000'0111u, 0b0u}, 9u);
    matrix.add_column({0b0111'0000u, 0b1u}, {0b1000'1111u, 0b0u}, 9u);
    matrix.add_column({0b0110'0000u}, {0b0000'1111u}, 7u);
    matrix.add_column({0b1110'0000u}, {0b0001'1111u}, 7u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        { -0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
        { -1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
        { -2, -1, -1, -1, -2, -3, -4, -5, -6, -7},
        { -3, -2, -2, -2, -2, -2, -3, -4, -5, -6},
        { -4, -3, -3, -3, -3, -3, -3, -3, -4, -5},
        { -5, -4, -3, -4, -4, -4, -4, -4, -4, -4},
        {INF, -5, -4, -3, -4, -5, -5, -5, -5, -5},
        {INF,INF, -5, -4, -4, -4, -5, -6,INF,INF},
        {INF,INF,INF, -5, -5, -5, -5, -5,INF,INF},
        {INF,INF,INF,INF,INF,INF,INF,INF,INF,INF}
    };

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, multiple_words_2)
{
    matrix_type<false, true> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u, 0b1111'1111u}, {0b0000'0000u, 0b0000'0000u}, 9u);
    matrix.add_column({0b1111'1110u, 0b1111'1111u}, {0b0000'0001u, 0b0000'0000u}, 10u);
    matrix.add_column({0b1111'1000u, 0b1111'1111u}, {0b0000'0001u, 0b0000'0000u}, 11u);
    matrix.add_column({0b1110'0000u, 0b1111'1111u}, {0b0000'0001u, 0b0000'0000u}, 12u);
    matrix.add_column({0b1000'0000u, 0b1111'1111u}, {0b0000'0001u, 0b0000'0000u}, 13u);
    matrix.add_column({0b0000'0100u, 0b1111'1110u}, {0b0000'0011u, 0b0000'0000u}, 14u);
    matrix.add_column({0b0001'1000u, 0b1111'1000u}, {0b0000'0111u, 0b0000'0000u}, 15u);
    matrix.add_column({0b0110'0000u, 0b1110'0000u}, {0b0000'0111u, 0b0000'0000u}, 16u);
    matrix.add_column({0b1000'0000u, 0b1000'0001u, 0b1u}, {0b0000'0111u, 0b0000'0000u, 0b0u}, 17u);
    matrix.add_column({0b0000'0000u, 0b0000'0110u, 0b0u}, {0b0000'0111u, 0b0000'0000u, 0b0u}, 18u);

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        { -0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
        { -1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
        { -2, -1, -1, -2, -3, -3, -4, -5, -6, -7},
        { -3, -2, -1, -2, -3, -4, -3, -4, -5, -6},
        { -4, -3, -2, -2, -3, -4, -4, -4, -5, -6},
        { -5, -4, -3, -2, -3, -4, -5, -4, -5, -6},
        { -6, -5, -4, -3, -3, -4, -5, -5, -5, -6},
        { -7, -6, -5, -4, -3, -4, -5, -6, -5, -6},
        { -8, -7, -6, -5, -4, -4, -5, -6, -6, -6},
        {INF, -8, -7, -6, -5, -4, -5, -6, -7, -6},
        {INF,INF, -8, -7, -6, -5, -5, -6, -7, -7},
        {INF,INF,INF, -8, -7, -6, -5, -6, -7, -8},
        {INF,INF,INF,INF, -8, -7, -6, -6, -7, -8},
        {INF,INF,INF,INF,INF, -8, -7, -6, -7, -8},
        {INF,INF,INF,INF,INF,INF, -8, -7, -7, -8},
        {INF,INF,INF,INF,INF,INF,INF, -8, -7, -8},
        {INF,INF,INF,INF,INF,INF,INF,INF, -8, -8},
        {INF,INF,INF,INF,INF,INF,INF,INF,INF, -8}
    };

    EXPECT_EQ(result, expect);
}
