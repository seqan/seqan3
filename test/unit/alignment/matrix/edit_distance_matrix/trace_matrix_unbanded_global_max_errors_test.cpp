// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/edit_distance_trace_matrix_full.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "edit_distance_trace_matrix.hpp"

TEST(global_max_errors, empty)
{
    matrix_type<false, true> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, epsilon)
{
    matrix_type<false, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, epsilon_row)
{
    matrix_type<false, true> matrix{1u};

    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 1u);
    matrix.add_column({}, {}, {}, 0u);
    matrix.add_column({}, {}, {}, 0u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, N, N}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_1)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 6u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 7u);
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u}, 8u);
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u}, 9u);
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u}, 9u);
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u}, 9u);
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u}, 9u);
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u}, 9u);
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u}, 7u);
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u}, 7u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l, l, l, l, l, l},
                                                                      {u, D, Dl, l, l, l, l, l, l, Dl},
                                                                      {u, u, D, D, Dl, l, l, l, l, l},
                                                                      {u, u, Du, Du, D, D, Dl, l, l, l},
                                                                      {u, u, Du, Du, Du, Du, D, D, Dl, l},
                                                                      {u, Du, D, Dul, Du, Du, Du, Du, D, D},
                                                                      {N, u, u, D, Dl, Dul, Du, Du, Du, Du},
                                                                      {N, N, u, u, D, D, Dl, Dul, N, N},
                                                                      {N, N, N, u, Du, Du, D, D, N, N}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_2)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 5u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 6u);
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u}, 7u);
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u}, 8u);
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u}, 8u);
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u}, 8u);
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u}, 6u);
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u}, 6u);
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u}, 6u);
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u}, 6u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l, l, l, l, l, l},
                                                                      {u, D, Dl, l, l, l, l, l, l, Dl},
                                                                      {u, u, D, D, Dl, l, l, l, l, l},
                                                                      {u, u, Du, Du, D, D, Dl, l, l, l},
                                                                      {u, u, Du, Du, Du, Du, D, D, Dl, l},
                                                                      {N, Du, D, Dul, Du, Du, Du, Du, D, D},
                                                                      {N, N, u, D, Dl, Dul, N, N, N, N},
                                                                      {N, N, N, u, D, D, N, N, N, N},
                                                                      {N, N, N, N, N, N, N, N, N, N}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, single_word_3)
{
    matrix_type<false, true> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u}, 4u);
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u}, 5u);
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u}, 6u);
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u}, 7u);
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u}, 5u);
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u}, 5u);
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u}, 5u);
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u}, 5u);
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u}, 0u);
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u}, 0u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l, l, l, l, N, N},
                                                                      {u, D, Dl, l, l, l, l, l, N, N},
                                                                      {u, u, D, D, Dl, l, l, l, N, N},
                                                                      {u, u, Du, Du, D, D, Dl, l, N, N},
                                                                      {N, u, Du, Du, Du, Du, D, D, N, N},
                                                                      {N, N, D, Dul, N, N, N, N, N, N},
                                                                      {N, N, N, D, N, N, N, N, N, N},
                                                                      {N, N, N, N, N, N, N, N, N, N},
                                                                      {N, N, N, N, N, N, N, N, N, N}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, multiple_words_1)
{
    matrix_type<false, true> matrix{10u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0u}, {0b0000'0000u, 0b0u}, {0b1111'1111u, 0b1u}, 6u);
    matrix.add_column({0b0000'0000u, 0b0u}, {0b0001'0001u, 0b1u}, {0b1111'1110u, 0b1u}, 7u);
    matrix.add_column({0b0000'0001u, 0b0u}, {0b0001'1111u, 0b1u}, {0b1110'1100u, 0b1u}, 8u);
    matrix.add_column({0b0001'0001u, 0b0u}, {0b0011'1110u, 0b0u}, {0b1101'1100u, 0b1u}, 9u);
    matrix.add_column({0b0010'0011u, 0b0u}, {0b1111'1110u, 0b1u}, {0b1001'1000u, 0b1u}, 9u);
    matrix.add_column({0b0010'0011u, 0b0u}, {0b1111'1100u, 0b1u}, {0b1011'1000u, 0b1u}, 9u);
    matrix.add_column({0b0100'0111u, 0b0u}, {0b1111'1100u, 0b1u}, {0b0011'0000u, 0b1u}, 9u);
    matrix.add_column({0b0100'0111u, 0b0u}, {0b1111'1000u, 0b1u}, {0b0111'0000u, 0b1u}, 9u);
    matrix.add_column({0b1000'1111u, 0b0u}, {0b1111'1000u, 0b1u}, {0b0110'0000u, 0b0u}, 7u);
    matrix.add_column({0b1000'1111u, 0b0u}, {0b1111'0001u, 0b1u}, {0b1110'0000u, 0b0u}, 7u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l, l, l, l, l, l},
                                                                      {u, D, Dl, l, l, l, l, l, l, Dl},
                                                                      {u, u, D, D, Dl, l, l, l, l, l},
                                                                      {u, u, Du, Du, D, D, Dl, l, l, l},
                                                                      {u, u, Du, Du, Du, Du, D, D, Dl, l},
                                                                      {u, Du, D, Dul, Du, Du, Du, Du, D, D},
                                                                      {N, u, u, D, Dl, Dul, Du, Du, Du, Du},
                                                                      {N, N, u, u, D, D, Dl, Dul, N, N},
                                                                      {N, N, N, u, Du, Du, D, D, N, N},
                                                                      {N, N, N, N, N, N, N, N, N, N}};

    EXPECT_EQ(result, expect);
}

TEST(global_max_errors, multiple_words_2)
{
    matrix_type<false, true> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1111'1111u, 0b1u},
                      9u);
    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0011u, 0b0000'0011u, 0b0u},
                      {0b1111'1110u, 0b1111'1111u, 0b1u},
                      10u);
    matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                      {0b0000'1110u, 0b0000'1100u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u},
                      11u);
    matrix.add_column({0b0000'0111u, 0b0000'0000u, 0b0u},
                      {0b0011'1110u, 0b0011'0000u, 0b0u},
                      {0b1110'0000u, 0b1111'1111u, 0b1u},
                      12u);
    matrix.add_column({0b0001'1111u, 0b0000'0000u, 0b0u},
                      {0b1111'1110u, 0b1100'0000u, 0b1u},
                      {0b1000'0000u, 0b1111'1111u, 0b1u},
                      13u);
    matrix.add_column({0b0111'1101u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b0000'0011u, 0b0u},
                      {0b0000'0100u, 0b1111'1110u, 0b1u},
                      14u);
    matrix.add_column({0b1111'0011u, 0b0000'0001u, 0b0u},
                      {0b1111'1100u, 0b0000'1111u, 0b0u},
                      {0b0001'1000u, 0b1111'1000u, 0b1u},
                      15u);
    matrix.add_column({0b1100'0111u, 0b0000'0111u, 0b0u},
                      {0b1111'1000u, 0b0011'1111u, 0b0u},
                      {0b0110'0000u, 0b1110'0000u, 0b1u},
                      16u);
    matrix.add_column({0b0001'1111u, 0b0001'1111u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u},
                      {0b1000'0000u, 0b1000'0001u, 0b1u},
                      17u);
    matrix.add_column({0b0111'1111u, 0b0111'1100u, 0b0u},
                      {0b1111'1011u, 0b1111'1111u, 0b1u},
                      {0b0000'0000u, 0b0000'0110u, 0b0u},
                      18u);

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l, l, l, l, l, l},
                                                                      {u, D, l, l, l, Dl, l, l, l, Dl},
                                                                      {u, Du, D, Dl, Dl, D, l, l, l, Dl},
                                                                      {u, u, D, Dl, Dl, Dul, D, l, l, l},
                                                                      {u, u, Du, D, Dl, Dl, Du, D, Dl, Dl},
                                                                      {u, u, u, D, Dl, Dl, Dul, D, Dl, Dl},
                                                                      {u, u, u, Du, D, Dl, Dl, Du, D, Dl},
                                                                      {u, u, u, u, D, Dl, Dl, Dul, D, Dl},
                                                                      {u, u, u, u, Du, D, Dl, Dl, Du, D},
                                                                      {N, Du, u, u, u, D, Dl, Dl, Dul, D},
                                                                      {N, N, u, u, u, Du, D, Dl, Dl, Du},
                                                                      {N, N, N, u, u, u, D, Dl, Dl, Dul},
                                                                      {N, N, N, N, u, u, Du, D, Dl, Dl},
                                                                      {N, N, N, N, N, u, u, D, Dl, Dl},
                                                                      {N, N, N, N, N, N, u, Du, D, Dl},
                                                                      {N, N, N, N, N, N, N, u, D, Dl},
                                                                      {N, N, N, N, N, N, N, N, Du, D},
                                                                      {N, N, N, N, N, N, N, N, N, D}};

    EXPECT_EQ(result, expect);
}
