// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/edit_distance_trace_matrix_full.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "edit_distance_trace_matrix.hpp"

TEST(global, empty)
{
    matrix_type<false, false> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(global, epsilon)
{
    matrix_type<false, false> matrix{1u};

    matrix.add_column({}, {}, {});

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N}};

    EXPECT_EQ(result, expect);
}

TEST(global, epsilon_row)
{
    matrix_type<false, false> matrix{1u};

    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});
    matrix.add_column({}, {}, {});

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l}};

    EXPECT_EQ(result, expect);
}

TEST(global, single_word)
{
    matrix_type<false, false> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u}, {0b0000'0000u}, {0b1111'1111u});
    matrix.add_column({0b0000'0000u}, {0b0001'0001u}, {0b1111'1110u});
    matrix.add_column({0b0000'0001u}, {0b0001'1111u}, {0b1110'1100u});
    matrix.add_column({0b0001'0001u}, {0b0011'1110u}, {0b1101'1100u});
    matrix.add_column({0b0010'0011u}, {0b1111'1110u}, {0b1001'1000u});
    matrix.add_column({0b0010'0011u}, {0b1111'1100u}, {0b1011'1000u});
    matrix.add_column({0b0100'0111u}, {0b1111'1100u}, {0b0011'0000u});
    matrix.add_column({0b0100'0111u}, {0b1111'1000u}, {0b0111'0000u});
    matrix.add_column({0b1000'1111u}, {0b1111'1000u}, {0b0110'0000u});
    matrix.add_column({0b1000'1111u}, {0b1111'0001u}, {0b1110'0000u});

    // row-wise matrix
    std::vector<std::vector<seqan3::detail::trace_directions>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<seqan3::detail::trace_directions>> expect{{N, l, l, l, l, l, l, l, l, l},
                                                                      {u, D, Dl, l, l, l, l, l, l, Dl},
                                                                      {u, u, D, D, Dl, l, l, l, l, l},
                                                                      {u, u, Du, Du, D, D, Dl, l, l, l},
                                                                      {u, u, Du, Du, Du, Du, D, D, Dl, l},
                                                                      {u, Du, D, Dul, Du, Du, Du, Du, D, D},
                                                                      {u, u, u, D, Dl, Dul, Du, Du, Du, Du},
                                                                      {u, u, u, u, D, D, Dl, Dul, Du, Du},
                                                                      {u, u, u, u, Du, Du, D, D, Dl, Dul}};

    EXPECT_EQ(result, expect);
}

TEST(global, multiple_words)
{
    matrix_type<false, false> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0000u, 0b0000'0000u, 0b0u},
                      {0b0000'0011u, 0b0000'0011u, 0b0u},
                      {0b1111'1110u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0001u, 0b0000'0000u, 0b0u},
                      {0b0000'1110u, 0b0000'1100u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0000'0111u, 0b0000'0000u, 0b0u},
                      {0b0011'1110u, 0b0011'0000u, 0b0u},
                      {0b1110'0000u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0001'1111u, 0b0000'0000u, 0b0u},
                      {0b1111'1110u, 0b1100'0000u, 0b1u},
                      {0b1000'0000u, 0b1111'1111u, 0b1u});
    matrix.add_column({0b0111'1101u, 0b0000'0000u, 0b0u},
                      {0b1111'1111u, 0b0000'0011u, 0b0u},
                      {0b0000'0100u, 0b1111'1110u, 0b1u});
    matrix.add_column({0b1111'0011u, 0b0000'0001u, 0b0u},
                      {0b1111'1100u, 0b0000'1111u, 0b0u},
                      {0b0001'1000u, 0b1111'1000u, 0b1u});
    matrix.add_column({0b1100'0111u, 0b0000'0111u, 0b0u},
                      {0b1111'1000u, 0b0011'1111u, 0b0u},
                      {0b0110'0000u, 0b1110'0000u, 0b1u});
    matrix.add_column({0b0001'1111u, 0b0001'1111u, 0b0u},
                      {0b1111'1000u, 0b1111'1111u, 0b1u},
                      {0b1000'0000u, 0b1000'0001u, 0b1u});
    matrix.add_column({0b0111'1111u, 0b0111'1100u, 0b0u},
                      {0b1111'1011u, 0b1111'1111u, 0b1u},
                      {0b0000'0000u, 0b0000'0110u, 0b0u});

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
                                                                      {u, Du, u, u, u, D, Dl, Dl, Dul, D},
                                                                      {u, Du, u, u, u, Du, D, Dl, Dl, Du},
                                                                      {u, u, Du, u, u, u, D, Dl, Dl, Dul},
                                                                      {u, u, Du, u, u, u, Du, D, Dl, Dl},
                                                                      {u, u, u, Du, u, u, u, D, Dl, Dl},
                                                                      {u, u, u, Du, u, u, u, Du, D, Dl},
                                                                      {u, u, u, u, Du, u, u, u, D, Dl},
                                                                      {u, u, u, u, Du, u, u, u, Du, D},
                                                                      {u, u, u, u, Du, u, u, u, Du, D}};

    EXPECT_EQ(result, expect);
}
