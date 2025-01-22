// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/detail/debug_matrix.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_same_type.hpp>

using seqan3::operator""_dna4;
using seqan3::operator|;

struct debug_matrix_test : public ::testing::Test
{
    static constexpr std::nullopt_t inf = std::nullopt;

    std::vector<seqan3::dna4> first_sequence = "AACACGTTAACCGGTT"_dna4;
    std::vector<seqan3::dna4> second_sequence = "ACGTACGT"_dna4;

    std::vector<bool> masking_matrix{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
                                     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0,
                                     1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
                                     0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<bool> transposed_masking_matrix{
        1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1,
        1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
        0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1,
        1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1,
        1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    std::vector<bool> masking_matrix_s9u_7u{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
                                            0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<bool> transposed_masking_matrix_s7u_9u{1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1,
                                                       1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
                                                       1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0};

    std::vector<int> scores{
        -0,  -1,  -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, -13, -14, -15, -16, -1,  -0,  -1, -2,  -3,
        -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, -13, -14, -15, -2,  -1,  -1,  -1,  -2,  -3,  -4,  -5, -6,  -7,
        -8,  -9,  -10, -11, -12, -13, -14, -3,  -2,  -2,  -2,  -2,  -3,  -3,  -4,  -5,  -6,  -7,  -8,  -9, -10, -11,
        -12, -13, -4,  -3,  -3,  -3,  -3,  -3,  -4,  -3,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, -5, -4,  -3,
        -4,  -3,  -4,  -4,  -4,  -4,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -6,  -5,  -4,  -3,  -4,  -3, -4,  -5,
        -5,  -5,  -5,  -5,  -6,  -7,  -8,  -9,  -10, -7,  -6,  -5,  -4,  -4,  -4,  -3,  -4,  -5,  -6,  -6, -6,  -6,
        -6,  -7,  -8,  -9,  -8,  -7,  -6,  -5,  -5,  -5,  -4,  -3,  -4,  -5,  -6,  -7,  -7,  -7,  -7,  -7, -8};

    seqan3::detail::row_wise_matrix<int> score_matrix{seqan3::detail::number_rows{9u},
                                                      seqan3::detail::number_cols{17u},
                                                      scores};

    seqan3::detail::row_wise_matrix<int> transposed_score_matrix{
        seqan3::detail::number_rows{17u},
        seqan3::detail::number_cols{9u},
        std::vector{-0,  -1, -2, -3, -4,  -5,  -6,  -7,  -8,  -1,  -0,  -1, -2, -3, -4, -5,  -6,  -7,  -2,  -1,
                    -1,  -2, -3, -3, -4,  -5,  -6,  -3,  -2,  -1,  -2,  -3, -4, -3, -4, -5,  -4,  -3,  -2,  -2,
                    -3,  -3, -4, -4, -5,  -5,  -4,  -3,  -3,  -3,  -4,  -3, -4, -5, -6, -5,  -4,  -3,  -4,  -4,
                    -4,  -3, -4, -7, -6,  -5,  -4,  -3,  -4,  -5,  -4,  -3, -8, -7, -6, -5,  -4,  -4,  -5,  -5,
                    -4,  -9, -8, -7, -6,  -5,  -4,  -5,  -6,  -5,  -10, -9, -8, -7, -6, -5,  -5,  -6,  -6,  -11,
                    -10, -9, -8, -7, -6,  -5,  -6,  -7,  -12, -11, -10, -9, -8, -7, -6, -6,  -7,  -13, -12, -11,
                    -10, -9, -8, -7, -6,  -7,  -14, -13, -12, -11, -10, -9, -8, -7, -7, -15, -14, -13, -12, -11,
                    -10, -9, -8, -7, -16, -15, -14, -13, -12, -11, -10, -9, -8}};

    seqan3::detail::row_wise_matrix<std::optional<int>> masked_score_matrix{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{17u},
        std::vector<std::optional<int>>{
            -0,  -1,  -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, -13, -14, -15, -16, -1,  -0,  -1,
            -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, -13, -14, -15, -2,  -1,  -1,  -1,  -2,  -3,
            -4,  -5,  -6,  -7,  -8,  -9,  -10, -11, -12, inf, -14, -3,  -2,  -2,  -2,  -2,  -3,  -3,  -4,  -5,
            -6,  -7,  -8,  -9,  -10, -11, inf, -13, -4,  -3,  -3,  -3,  -3,  -3,  -4,  -3,  -4,  -5,  -6,  -7,
            -8,  -9,  -10, inf, -12, inf, -4,  -3,  -4,  -3,  -4,  -4,  inf, inf, inf, -5,  -6,  -7,  -8,  -9,
            inf, -11, inf, inf, -4,  -3,  -4,  -3,  inf, inf, inf, inf, inf, -5,  -6,  -7,  inf, inf, -10, inf,
            inf, inf, -4,  -4,  inf, inf, inf, inf, inf, inf, inf, -6,  -6,  inf, inf, -9,  inf, inf, inf, inf,
            inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, inf, -8}};

    seqan3::detail::row_wise_matrix<int> score_matrix_s9u_7u{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{7u},
        std::vector{-0, -1, -2, -3, -4, -5, -6, -1, -0, -1, -2, -3, -4, -5, -2, -1, -1, -1, -2, -3, -4,
                    -3, -2, -2, -2, -2, -3, -3, -4, -3, -3, -3, -3, -3, -4, -5, -4, -3, -4, -3, -4, -4,
                    -6, -5, -4, -3, -4, -3, -4, -7, -6, -5, -4, -4, -4, -3, -8, -7, -6, -5, -5, -5, -4}};

    seqan3::detail::row_wise_matrix<int> transposed_score_matrix_s9u_7u{
        seqan3::detail::number_rows{7u},
        seqan3::detail::number_cols{9u},
        std::vector{-0, -1, -2, -3, -4, -5, -6, -7, -8, -1, -0, -1, -2, -3, -4, -5, -6, -7, -2, -1, -1,
                    -2, -3, -3, -4, -5, -6, -3, -2, -1, -2, -3, -4, -3, -4, -5, -4, -3, -2, -2, -3, -3,
                    -4, -4, -5, -5, -4, -3, -3, -3, -4, -3, -4, -5, -6, -5, -4, -3, -4, -4, -4, -3, -4}};

    seqan3::detail::row_wise_matrix<int> score_matrix_s4u_17u{
        seqan3::detail::number_rows{4u},
        seqan3::detail::number_cols{17u},
        std::vector{-0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16,
                    -1, -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,  -10, -11, -12, -13, -14, -15,
                    -2, -1, -1, -1, -2, -3, -4, -5, -6, -7, -8,  -9,  -10, -11, -12, -13, -14,
                    -3, -2, -2, -2, -2, -3, -3, -4, -5, -6, -7,  -8,  -9,  -10, -11, -12, -13}};

    seqan3::detail::trace_directions N{}, D{seqan3::detail::trace_directions::diagonal},
        L{seqan3::detail::trace_directions::left}, U{seqan3::detail::trace_directions::up}, DL{D | L}, DU{D | U},
        UL{U | L}, DUL{D | U | L};

    std::vector<seqan3::detail::trace_directions> traces{
        N,  L,  L,  L,  L,  L,   L,  L,  L,  L, L,   L,  L,   L,  L,  L, L, U,   D,  DL, L,  DL, L,   L, L,   L,
        DL, DL, L,  L,  L,  L,   L,  L,  U,  U, D,   D,  L,   DL, L,  L, L, L,   L,  DL, DL, L,  L,   L, L,   U,
        U,  DU, DU, D,  DL, D,   L,  L,  L,  L, L,   L,  DL,  DL, L,  L, U, U,   DU, DU, DU, D,  DUL, D, DL,  L,
        L,  L,  L,  L,  L,  DL,  DL, U,  DU, D, DUL, D,  DUL, D,  U,  D, D, DL,  L,  L,  L,  L,  L,   L, U,   U,
        U,  D,  UL, D,  L,  DUL, DU, DU, D,  D, DL,  L,  L,   L,  L,  U, U, U,   U,  D,  U,  D,  L,   L, DUL, DU,
        DU, D,  D,  DL, L,  L,   U,  U,  U,  U, DU,  DU, U,   D,  DL, L, L, DUL, DU, DU, D,  D,  DL};

    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> trace_matrix{seqan3::detail::number_rows{9u},
                                                                                   seqan3::detail::number_cols{17u},
                                                                                   traces};

    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> transposed_trace_matrix{
        seqan3::detail::number_rows{17u},
        seqan3::detail::number_cols{9u},
        std::vector{N,   L,  L, L,   L,  L, L,   L,  L,  U,  D,   L, L,  L,  DL,  L,  L,  L, U,  DU,  D,  DL,
                    DL,  D,  L, L,   L,  U, U,   D,  DL, DL, DUL, D, L,  L,  U,   DU, U,  D, DL, D,   UL, D,
                    DL,  U,  U, DU,  DU, D, DUL, D,  L,  DL, U,   U, U,  D,  DUL, D,  U,  D, L,  U,   U,  U,
                    U,   D,  L, DUL, U,  D, U,   U,  U,  U,  DU,  D, DL, U,  DU,  U,  DU, U, U,  U,   D,  DL,
                    DUL, U,  U, DU,  U,  U, U,   DU, D,  DL, U,   U, U,  DU, U,   U,  U,  D, DL, DUL, U,  U,
                    DU,  U,  U, U,   DU, D, DL,  U,  U,  U,  DU,  U, U,  U,  D,   DL, U,  U, U,  DU,  U,  U,
                    U,   DU, D, U,   U,  U, U,   DU, U,  U,  U,   D, U,  U,  U,   U,  DU, U, U,  U,   DU}};

    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> masked_trace_matrix{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{17u},
        std::vector{N,  L,  L,  L, L,  L, L,  L, L,  L, L,   L, L,   L,  L, L, L, U,  D,  DL, L,  DL, L,   L, L,  L,
                    DL, DL, L,  L, L,  L, L,  L, U,  U, D,   D, L,   DL, L, L, L, L,  L,  DL, DL, L,  L,   N, L,  U,
                    U,  DU, DU, D, DL, D, L,  L, L,  L, L,   L, DL,  DL, N, L, U, U,  DU, DU, DU, D,  DUL, D, DL, L,
                    L,  L,  L,  L, L,  N, DL, N, DU, D, DUL, D, DUL, D,  N, N, N, DL, L,  L,  L,  L,  N,   L, N,  N,
                    U,  D,  UL, D, N,  N, N,  N, N,  D, DL,  L, N,   N,  L, N, N, N,  U,  D,  N,  N,  N,   N, N,  N,
                    N,  D,  D,  N, N,  L, N,  N, N,  N, N,   N, N,   N,  N, N, N, N,  N,  N,  N,  N,  DL}};

    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> trace_matrix_s9u_7u{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{7u},
        std::vector{N, L, L,  L,  L,  L,  L, U, D, DL, L,  DL, L, L,   U, U,  D, D,   L,  DL,  L,
                    U, U, DU, DU, D,  DL, D, U, U, DU, DU, DU, D, DUL, U, DU, D, DUL, D,  DUL, D,
                    U, U, U,  D,  UL, D,  L, U, U, U,  U,  D,  U, D,   U, U,  U, U,   DU, DU,  U}};

    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> trace_matrix_s4u_17u{
        seqan3::detail::number_rows{4u},
        seqan3::detail::number_cols{17u},
        std::vector{N,  L, L, L,  L,  L, L, L,  L,  L, L,  L, L, L, L, L, L,  U, D,  DL, L, DL, L,
                    L,  L, L, DL, DL, L, L, L,  L,  L, L,  U, U, D, D, L, DL, L, L,  L,  L, L,  DL,
                    DL, L, L, L,  L,  U, U, DU, DU, D, DL, D, L, L, L, L, L,  L, DL, DL, L, L}};

    template <typename score_matrix_t>
    void test_score_matrix(score_matrix_t && matrix)
    {
        EXPECT_EQ(matrix.cols(), 17u);
        EXPECT_EQ(matrix.rows(), 9u);

        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{0u}}), -0);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{6u}}), -6);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{16u}}), -16);

        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{3u}, seqan3::detail::column_index_type{0u}}), -3);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{3u}, seqan3::detail::column_index_type{6u}}), -3);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{3u}, seqan3::detail::column_index_type{16u}}), -13);

        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{4u}, seqan3::detail::column_index_type{0u}}), -4);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{4u}, seqan3::detail::column_index_type{6u}}), -4);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{4u}, seqan3::detail::column_index_type{16u}}), -12);

        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{8u}, seqan3::detail::column_index_type{0u}}), -8);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{8u}, seqan3::detail::column_index_type{6u}}), -4);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{8u}, seqan3::detail::column_index_type{16u}}), -8);

        for (size_t row = 0; row < matrix.rows(); row++)
            for (size_t col = 0; col < matrix.cols(); col++)
                EXPECT_EQ((matrix.at({seqan3::detail::row_index_type{row}, seqan3::detail::column_index_type{col}})),
                          scores[row * matrix.cols() + col]);
    }

    template <typename trace_matrix_t>
    void test_trace_matrix(trace_matrix_t && matrix)
    {
        EXPECT_EQ(matrix.cols(), 17u);
        EXPECT_EQ(matrix.rows(), 9u);

        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{0u}}), N);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{3u}, seqan3::detail::column_index_type{6u}}), D);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{3u}, seqan3::detail::column_index_type{0u}}), U);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{6u}}), L);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{8u}, seqan3::detail::column_index_type{5u}}), DU);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{5u}}), DL);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{6u}, seqan3::detail::column_index_type{4u}}), UL);
        EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{4u}, seqan3::detail::column_index_type{6u}}), DUL);

        for (size_t row = 0; row < matrix.rows(); row++)
            for (size_t col = 0; col < matrix.cols(); col++)
                EXPECT_EQ(matrix.at({seqan3::detail::row_index_type{row}, seqan3::detail::column_index_type{col}}),
                          traces[row * matrix.cols() + col]);
    }
};

using score_matrix_test = debug_matrix_test;
using trace_matrix_test = debug_matrix_test;

template <typename>
struct debug_matrix_traits;

template <seqan3::detail::matrix matrix_t, typename first_sequence_t, typename second_sequence_t>
struct debug_matrix_traits<seqan3::detail::debug_matrix<matrix_t, first_sequence_t, second_sequence_t>>
{
    using matrix_type = matrix_t;
    using first_sequence_type = first_sequence_t;
    using second_sequence_type = second_sequence_t;
};

#pragma GCC diagnostic push
// Ignore bogus warnings in fortified gcc13 build
// _FORTIFY_SOURCE=2 may cause conforming programs to fail
#if not defined(__clang__) && defined(_FORTIFY_SOURCE) && (_FORTIFY_SOURCE == 2)
#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif

TEST_F(debug_matrix_test, matrix_concept)
{
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::row_wise_matrix<int>>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::row_wise_matrix<int> &>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::row_wise_matrix<int> const>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::row_wise_matrix<int> const &>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::debug_matrix<seqan3::detail::row_wise_matrix<int>>>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::debug_matrix<seqan3::detail::row_wise_matrix<int> &>>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::debug_matrix<seqan3::detail::row_wise_matrix<int> const>>));
    EXPECT_TRUE((seqan3::detail::matrix<seqan3::detail::debug_matrix<seqan3::detail::row_wise_matrix<int> const &>>));
}

TEST_F(debug_matrix_test, construct_with_references)
{
    using debug_matrix_type = decltype(seqan3::detail::debug_matrix{score_matrix, first_sequence, second_sequence});
    using matrix_type = typename debug_matrix_traits<debug_matrix_type>::matrix_type;
    using first_sequence_type = typename debug_matrix_traits<debug_matrix_type>::first_sequence_type;
    using second_sequence_type = typename debug_matrix_traits<debug_matrix_type>::second_sequence_type;

    EXPECT_SAME_TYPE(matrix_type, seqan3::detail::row_wise_matrix<int> &);
    EXPECT_SAME_TYPE(first_sequence_type, std::vector<seqan3::dna4> &);
    EXPECT_SAME_TYPE(second_sequence_type, std::vector<seqan3::dna4> &);
}

TEST_F(debug_matrix_test, construct_with_move)
{
    using debug_matrix_type = decltype(seqan3::detail::debug_matrix{std::move(score_matrix),
                                                                    std::move(first_sequence),
                                                                    std::move(second_sequence)});
    using matrix_type = typename debug_matrix_traits<debug_matrix_type>::matrix_type;
    using first_sequence_type = typename debug_matrix_traits<debug_matrix_type>::first_sequence_type;
    using second_sequence_type = typename debug_matrix_traits<debug_matrix_type>::second_sequence_type;

    EXPECT_SAME_TYPE(matrix_type, seqan3::detail::row_wise_matrix<int>);
    EXPECT_SAME_TYPE(first_sequence_type, std::vector<seqan3::dna4>);
    EXPECT_SAME_TYPE(second_sequence_type, std::vector<seqan3::dna4>);
}

TEST_F(score_matrix_test, other_matrix)
{
    seqan3::detail::debug_matrix matrix{score_matrix};

    test_score_matrix(std::move(matrix));
}

TEST_F(score_matrix_test, sequences_other_matrix)
{
    seqan3::detail::debug_matrix matrix{score_matrix, first_sequence, second_sequence};

    EXPECT_EQ(matrix.first_sequence(), first_sequence);
    EXPECT_EQ(matrix.second_sequence(), second_sequence);

    test_score_matrix(std::move(matrix));
}

TEST_F(score_matrix_test, equal)
{
    // last entry of second row
    std::vector<int> scores_unequal{scores};
    scores_unequal[2 * 16] = -16;
    seqan3::detail::row_wise_matrix<int> score_matrix_unequal{seqan3::detail::number_rows{9u},
                                                              seqan3::detail::number_cols{17u},
                                                              std::move(scores_unequal)};

    seqan3::detail::debug_matrix matrix{score_matrix};

    EXPECT_EQ(matrix, score_matrix);
    EXPECT_EQ(matrix, matrix);
    EXPECT_FALSE(matrix == score_matrix_s9u_7u);
    EXPECT_FALSE(matrix == score_matrix_s4u_17u);
    EXPECT_FALSE(matrix == score_matrix_unequal);
}

TEST_F(score_matrix_test, not_equal)
{
    // last entry of second row
    std::vector<int> scores_unequal{scores};
    scores_unequal[2 * 16] = -16;
    seqan3::detail::row_wise_matrix<int> score_matrix_unequal{seqan3::detail::number_rows{9u},
                                                              seqan3::detail::number_cols{17u},
                                                              std::move(scores_unequal)};

    seqan3::detail::debug_matrix matrix{score_matrix};

    EXPECT_FALSE(matrix != score_matrix);
    EXPECT_FALSE(matrix != matrix);
    EXPECT_NE(matrix, score_matrix_s9u_7u);
    EXPECT_NE(matrix, score_matrix_s4u_17u);
    EXPECT_NE(matrix, score_matrix_unequal);
}

TEST_F(score_matrix_test, sub_matrix_lvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{score_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix sub_matrix = matrix.sub_matrix(9u, 7u);

    EXPECT_EQ(sub_matrix.rows(), 9u);
    EXPECT_EQ(sub_matrix.cols(), 7u);
    EXPECT_EQ(sub_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(sub_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(sub_matrix, score_matrix_s9u_7u);
}

TEST_F(score_matrix_test, sub_matrix_rvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{score_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix sub_matrix = std::move(matrix).sub_matrix(9u, 7u);

    EXPECT_EQ(sub_matrix.rows(), 9u);
    EXPECT_EQ(sub_matrix.cols(), 7u);
    EXPECT_EQ(sub_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(sub_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(sub_matrix, score_matrix_s9u_7u);
    EXPECT_EQ((seqan3::detail::debug_matrix{score_matrix}.sub_matrix(4u, 17u)), score_matrix_s4u_17u);
}

TEST_F(score_matrix_test, mask_matrix_lvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{score_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix mask_matrix = matrix.mask_matrix(masking_matrix);

    EXPECT_EQ(mask_matrix.rows(), 9u);
    EXPECT_EQ(mask_matrix.cols(), 17u);
    EXPECT_EQ(mask_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(mask_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(mask_matrix, masked_score_matrix);
}

TEST_F(score_matrix_test, mask_matrix_rvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{score_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix mask_matrix = std::move(matrix).mask_matrix(masking_matrix);

    EXPECT_EQ(mask_matrix.rows(), 9u);
    EXPECT_EQ(mask_matrix.cols(), 17u);
    EXPECT_EQ(mask_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(mask_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(mask_matrix, masked_score_matrix);
}

TEST_F(score_matrix_test, transpose_matrix_lvalue)
{
    auto first_sequence_expect = second_sequence;
    auto second_sequence_expect = first_sequence;
    seqan3::detail::debug_matrix matrix{score_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix transpose_matrix = matrix.transpose_matrix();

    EXPECT_EQ(transpose_matrix.rows(), 17u);
    EXPECT_EQ(transpose_matrix.cols(), 9u);
    EXPECT_EQ(transpose_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(transpose_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(transpose_matrix, transposed_score_matrix);
}

TEST_F(score_matrix_test, transpose_matrix_rvalue)
{
    auto first_sequence_expect = second_sequence;
    auto second_sequence_expect = first_sequence;
    seqan3::detail::debug_matrix matrix{score_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix transpose_matrix = std::move(matrix).transpose_matrix();

    EXPECT_EQ(transpose_matrix.rows(), 17u);
    EXPECT_EQ(transpose_matrix.cols(), 9u);
    EXPECT_EQ(transpose_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(transpose_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(transpose_matrix, transposed_score_matrix);
}

TEST_F(score_matrix_test, combine_sub_transpose_operations)
{
    seqan3::detail::debug_matrix matrix1{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix2{score_matrix, first_sequence, second_sequence};
    matrix1.transpose_matrix().sub_matrix(7u, 9u);
    matrix2.sub_matrix(9u, 7u).transpose_matrix();

    EXPECT_EQ(matrix1.rows(), 7u);
    EXPECT_EQ(matrix1.cols(), 9u);
    EXPECT_EQ(matrix2.rows(), 7u);
    EXPECT_EQ(matrix2.cols(), 9u);
    EXPECT_EQ(matrix1, transposed_score_matrix_s9u_7u);
    EXPECT_EQ(matrix2, transposed_score_matrix_s9u_7u);
    EXPECT_EQ(matrix1, matrix2);
}

TEST_F(score_matrix_test, combine_mask_transpose_operations)
{
    seqan3::detail::debug_matrix matrix1{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix2{score_matrix, first_sequence, second_sequence};
    matrix1.mask_matrix(masking_matrix).sub_matrix(9u, 7u);
    matrix2.sub_matrix(9u, 7u).mask_matrix(masking_matrix_s9u_7u);

    EXPECT_EQ(matrix1.rows(), 9u);
    EXPECT_EQ(matrix1.cols(), 7u);
    EXPECT_EQ(matrix2.rows(), 9u);
    EXPECT_EQ(matrix2.cols(), 7u);
    EXPECT_EQ(matrix1, matrix2);
}

TEST_F(score_matrix_test, combine_sub_mask_transpose_operations)
{
    seqan3::detail::debug_matrix matrix1{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix2{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix3{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix4{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix5{score_matrix, first_sequence, second_sequence};
    seqan3::detail::debug_matrix matrix6{score_matrix, first_sequence, second_sequence};
    matrix1.mask_matrix(masking_matrix).transpose_matrix().sub_matrix(7u, 9u);
    matrix2.mask_matrix(masking_matrix).sub_matrix(9u, 7u).transpose_matrix();
    matrix3.sub_matrix(9u, 7u).mask_matrix(masking_matrix_s9u_7u).transpose_matrix();
    matrix4.sub_matrix(9u, 7u).transpose_matrix().mask_matrix(transposed_masking_matrix_s7u_9u);
    matrix5.transpose_matrix().sub_matrix(7u, 9u).mask_matrix(transposed_masking_matrix_s7u_9u);
    matrix6.transpose_matrix().mask_matrix(transposed_masking_matrix).sub_matrix(7u, 9u);

    EXPECT_EQ(matrix1.rows(), 7u);
    EXPECT_EQ(matrix1.cols(), 9u);
    EXPECT_EQ(matrix2.rows(), 7u);
    EXPECT_EQ(matrix2.cols(), 9u);
    EXPECT_EQ(matrix3.rows(), 7u);
    EXPECT_EQ(matrix3.cols(), 9u);
    EXPECT_EQ(matrix4.rows(), 7u);
    EXPECT_EQ(matrix4.cols(), 9u);
    EXPECT_EQ(matrix5.rows(), 7u);
    EXPECT_EQ(matrix5.cols(), 9u);
    EXPECT_EQ(matrix6.rows(), 7u);
    EXPECT_EQ(matrix6.cols(), 9u);
    EXPECT_EQ(matrix1, matrix2);
    EXPECT_EQ(matrix1, matrix3);
    EXPECT_EQ(matrix1, matrix4);
    EXPECT_EQ(matrix1, matrix5);
    EXPECT_EQ(matrix1, matrix6);
    EXPECT_EQ(matrix2, matrix3);
    EXPECT_EQ(matrix2, matrix4);
    EXPECT_EQ(matrix2, matrix5);
    EXPECT_EQ(matrix2, matrix6);
    EXPECT_EQ(matrix3, matrix4);
    EXPECT_EQ(matrix3, matrix5);
    EXPECT_EQ(matrix3, matrix6);
    EXPECT_EQ(matrix4, matrix5);
    EXPECT_EQ(matrix4, matrix6);
    EXPECT_EQ(matrix5, matrix6);
}

TEST_F(trace_matrix_test, other_matrix)
{
    seqan3::detail::debug_matrix matrix{trace_matrix};

    test_trace_matrix(std::move(matrix));
}

TEST_F(trace_matrix_test, sequences_other_matrix)
{
    seqan3::detail::debug_matrix matrix{trace_matrix, first_sequence, second_sequence};

    EXPECT_EQ(matrix.first_sequence(), first_sequence);
    EXPECT_EQ(matrix.second_sequence(), second_sequence);

    test_trace_matrix(std::move(matrix));
}

TEST_F(trace_matrix_test, equal)
{
    // last entry of second row
    std::vector<seqan3::detail::trace_directions> traces_unequal{traces};
    traces_unequal[2 * 16] = seqan3::detail::trace_directions::up;
    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> trace_matrix_unequal{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{17u},
        std::move(traces_unequal)};

    seqan3::detail::debug_matrix matrix{trace_matrix};

    EXPECT_EQ(matrix, trace_matrix);
    EXPECT_EQ(matrix, matrix);
    EXPECT_FALSE(matrix == trace_matrix_s9u_7u);
    EXPECT_FALSE(matrix == trace_matrix_s4u_17u);
    EXPECT_FALSE(matrix == trace_matrix_unequal);
}

TEST_F(trace_matrix_test, not_equal)
{
    // last entry of second row
    std::vector<seqan3::detail::trace_directions> traces_unequal{traces};
    traces_unequal[2 * 16] = seqan3::detail::trace_directions::up;
    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> trace_matrix_unequal{
        seqan3::detail::number_rows{9u},
        seqan3::detail::number_cols{17u},
        std::move(traces_unequal)};

    seqan3::detail::debug_matrix matrix{trace_matrix};

    EXPECT_FALSE(matrix != trace_matrix);
    EXPECT_FALSE(matrix != matrix);
    EXPECT_NE(matrix, trace_matrix_s9u_7u);
    EXPECT_NE(matrix, trace_matrix_s4u_17u);
    EXPECT_NE(matrix, trace_matrix_unequal);
}

TEST_F(trace_matrix_test, sub_matrix_lvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{trace_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix sub_matrix = matrix.sub_matrix(9u, 7u);

    EXPECT_EQ(sub_matrix.rows(), 9u);
    EXPECT_EQ(sub_matrix.cols(), 7u);
    EXPECT_EQ(sub_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(sub_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(sub_matrix, trace_matrix_s9u_7u);
}

TEST_F(trace_matrix_test, sub_matrix_rvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{trace_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix sub_matrix = std::move(matrix).sub_matrix(9u, 7u);

    EXPECT_EQ(sub_matrix.rows(), 9u);
    EXPECT_EQ(sub_matrix.cols(), 7u);
    EXPECT_EQ(sub_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(sub_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(sub_matrix, trace_matrix_s9u_7u);
    EXPECT_EQ((seqan3::detail::debug_matrix{trace_matrix}.sub_matrix(4u, 17u)), trace_matrix_s4u_17u);
}

TEST_F(trace_matrix_test, mask_matrix_lvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{trace_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix mask_matrix = matrix.mask_matrix(masking_matrix);

    EXPECT_EQ(mask_matrix.rows(), 9u);
    EXPECT_EQ(mask_matrix.cols(), 17u);
    EXPECT_EQ(mask_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(mask_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(mask_matrix, masked_trace_matrix);
}

TEST_F(trace_matrix_test, mask_matrix_rvalue)
{
    auto first_sequence_expect = first_sequence;
    auto second_sequence_expect = second_sequence;
    seqan3::detail::debug_matrix matrix{trace_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix mask_matrix = std::move(matrix).mask_matrix(masking_matrix);

    EXPECT_EQ(mask_matrix.rows(), 9u);
    EXPECT_EQ(mask_matrix.cols(), 17u);
    EXPECT_EQ(mask_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(mask_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(mask_matrix, masked_trace_matrix);
}

TEST_F(trace_matrix_test, transpose_matrix_lvalue)
{
    auto first_sequence_expect = second_sequence;
    auto second_sequence_expect = first_sequence;
    seqan3::detail::debug_matrix matrix{trace_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix transpose_matrix = matrix.transpose_matrix();

    EXPECT_EQ(transpose_matrix.rows(), 17u);
    EXPECT_EQ(transpose_matrix.cols(), 9u);
    EXPECT_EQ(transpose_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(transpose_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(transpose_matrix, transposed_trace_matrix);
}

TEST_F(trace_matrix_test, transpose_matrix_rvalue)
{
    auto first_sequence_expect = second_sequence;
    auto second_sequence_expect = first_sequence;
    seqan3::detail::debug_matrix matrix{trace_matrix, std::move(first_sequence), std::move(second_sequence)};
    seqan3::detail::debug_matrix transpose_matrix = std::move(matrix).transpose_matrix();

    EXPECT_EQ(transpose_matrix.rows(), 17u);
    EXPECT_EQ(transpose_matrix.cols(), 9u);
    EXPECT_EQ(transpose_matrix.first_sequence(), first_sequence_expect);
    EXPECT_EQ(transpose_matrix.second_sequence(), second_sequence_expect);

    EXPECT_EQ(transpose_matrix, transposed_trace_matrix);
}

#pragma GCC diagnostic pop
