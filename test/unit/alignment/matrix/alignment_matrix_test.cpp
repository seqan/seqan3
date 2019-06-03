// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;
using namespace seqan3::detail;

struct no_config
{};

struct score_matrix : public ::testing::Test
{
    std::vector<dna4> database = "AACACGTTAACCGGTT"_dna4;
    std::vector<dna4> query = "ACGTACGT"_dna4;
    std::vector<int> scores
    {
       -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,
       -1, -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,
       -2, -1, -1, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,
       -3, -2, -2, -2, -2, -3, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,
       -4, -3, -3, -3, -3, -3, -4, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,
       -5, -4, -3, -4, -3, -4, -4, -4, -4, -4, -5, -6, -7, -8, -9,-10,-11,
       -6, -5, -4, -3, -4, -3, -4, -5, -5, -5, -5, -5, -6, -7, -8, -9,-10,
       -7, -6, -5, -4, -4, -4, -3, -4, -5, -6, -6, -6, -6, -6, -7, -8, -9,
       -8, -7, -6, -5, -5, -5, -4, -3, -4, -5, -6, -7, -7, -7, -7, -7, -8
    };
    std::vector<int> scores_shorter_cols
    {
       -0, -1, -2, -3, -4, -5, -6,
       -1, -0, -1, -2, -3, -4, -5,
       -2, -1, -1, -1, -2, -3, -4,
       -3, -2, -2, -2, -2, -3, -3,
       -4, -3, -3, -3, -3, -3, -4,
       -5, -4, -3, -4, -3, -4, -4,
       -6, -5, -4, -3, -4, -3, -4,
       -7, -6, -5, -4, -4, -4, -3,
       -8, -7, -6, -5, -5, -5, -4
    };
    std::vector<int> scores_shorter_rows
    {
       -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,
       -1, -0, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,-15,
       -2, -1, -1, -1, -2, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13,-14,
       -3, -2, -2, -2, -2, -3, -3, -4, -5, -6, -7, -8, -9,-10,-11,-12,-13
    };

    trace_directions N{},
        D{trace_directions::diagonal},
        L{trace_directions::left},
        U{trace_directions::up},
        DL{D|L}, DU{D|U}, UL{U|L}, DUL{D|U|L};

    std::vector<trace_directions> traces
    {
        N,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,
        U,  D,  DL, L,  DL, L,  L,  L,  L,  DL, DL, L,  L,  L,  L,  L,  L,
        U,  U,  D,  D,  L,  DL, L,  L,  L,  L,  L,  DL, DL, L,  L,  L,  L,
        U,  U,  DU, DU, D,  DL, D,  L,  L,  L,  L,  L,  L,  DL, DL, L,  L,
        U,  U,  DU, DU, DU, D,  DUL,D,  DL, L,  L,  L,  L,  L,  L,  DL, DL,
        U,  DU, D,  DUL,D,  DUL,D,  U,  D,  D,  DL, L,  L,  L,  L,  L,  L,
        U,  U,  U,  D,  UL, D,  L,  DUL,DU, DU, D,  D,  DL, L,  L,  L,  L,
        U,  U,  U,  U,  D,  U,  D,  L,  L,  DUL,DU, DU, D,  D,  DL, L,  L,
        U,  U,  U,  U,  DU, DU, U,  D,  DL, L,  L,  DUL,DU, DU, D,  D,  DL
    };
    std::vector<trace_directions> traces_shorter_rows
    {
        N,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,
        U,  D,  DL, L,  DL, L,  L,  L,  L,  DL, DL, L,  L,  L,  L,  L,  L,
        U,  U,  D,  D,  L,  DL, L,  L,  L,  L,  L,  DL, DL, L,  L,  L,  L,
        U,  U,  DU, DU, D,  DL, D,  L,  L,  L,  L,  L,  L,  DL, DL, L,  L
    };
    std::vector<trace_directions> traces_shorter_cols
    {
        N,  L,  L,  L,  L,  L,  L,
        U,  D,  DL, L,  DL, L,  L,
        U,  U,  D,  D,  L,  DL, L,
        U,  U,  DU, DU, D,  DL, D,
        U,  U,  DU, DU, DU, D,  DUL,
        U,  DU, D,  DUL,D,  DUL,D,
        U,  U,  U,  D,  UL, D,  L,
        U,  U,  U,  U,  D,  U,  D,
        U,  U,  U,  U,  DU, DU, U
    };

    template <typename ...options_t>
    void score_matrix_test(alignment_score_matrix<options_t...> matrix)
    {
        EXPECT_EQ(matrix.cols(), 17u);
        EXPECT_EQ(matrix.rows(), 9u);

        EXPECT_EQ(matrix.at(0, 0), -0);
        EXPECT_EQ(matrix.at(0, 6), -6);
        EXPECT_EQ(matrix.at(0, 16), -16);

        EXPECT_EQ(matrix.at(3, 0), -3);
        EXPECT_EQ(matrix.at(3, 6), -3);
        EXPECT_EQ(matrix.at(3, 16), -13);

        EXPECT_EQ(matrix.at(4, 0), -4);
        EXPECT_EQ(matrix.at(4, 6), -4);
        EXPECT_EQ(matrix.at(4, 16), -12);

        EXPECT_EQ(matrix.at(8, 0), -8);
        EXPECT_EQ(matrix.at(8, 6), -4);
        EXPECT_EQ(matrix.at(8, 16), -8);

        for (size_t row = 0; row < matrix.rows(); row++)
            for (size_t col = 0; col < matrix.cols(); col++)
                EXPECT_EQ(matrix.at(row, col), scores[row * matrix.cols() + col]);
    }

    template <typename ...options_t>
    void trace_matrix_test(alignment_trace_matrix<options_t...> matrix)
    {
        EXPECT_EQ(matrix.cols(), 17u);
        EXPECT_EQ(matrix.rows(), 9u);

        EXPECT_EQ(matrix.at(0, 0), N);
        EXPECT_EQ(matrix.at(3, 6), D);
        EXPECT_EQ(matrix.at(3, 0), U);
        EXPECT_EQ(matrix.at(0, 6), L);
        EXPECT_EQ(matrix.at(8, 5), DU);
        EXPECT_EQ(matrix.at(2, 5), DL);
        EXPECT_EQ(matrix.at(6, 4), UL);
        EXPECT_EQ(matrix.at(4, 6), DUL);

        for (size_t row = 0; row < matrix.rows(); row++)
            for (size_t col = 0; col < matrix.cols(); col++)
                EXPECT_EQ(matrix.at(row, col), traces[row * matrix.cols() + col]);
    }
};

using trace_matrix = score_matrix;

TEST_F(score_matrix, vector)
{
    alignment_score_matrix matrix{scores, 9u, 17u};

    score_matrix_test(std::move(matrix));
}

TEST_F(score_matrix, equal)
{
    // last entry of second row
    std::vector<int> scores_unequal{scores};
    scores_unequal[2 * 16] = -16;

    alignment_score_matrix matrix_shorter_cols{scores_shorter_cols, 9u, 7u};
    alignment_score_matrix matrix_shorter_rows{scores_shorter_rows, 4u, 17u};
    alignment_score_matrix matrix_unequal{std::move(scores_unequal), 9u, 17u};

    alignment_score_matrix matrix_vector{scores, 9u, 17u};

    EXPECT_EQ(matrix_vector, matrix_vector);
    EXPECT_FALSE(matrix_vector == matrix_shorter_cols);
    EXPECT_FALSE(matrix_vector == matrix_shorter_rows);
    EXPECT_FALSE(matrix_vector == matrix_unequal);
}

TEST_F(score_matrix, not_equal)
{
    // last entry of second row
    std::vector<int> scores_unequal{scores};
    scores_unequal[2 * 16] = -16;

    alignment_score_matrix matrix_shorter_cols{scores_shorter_cols, 9u, 7u};
    alignment_score_matrix matrix_shorter_rows{scores_shorter_rows, 4u, 17u};
    alignment_score_matrix matrix_unequal{std::move(scores_unequal), 9u, 17u};

    alignment_score_matrix matrix_vector{scores, 9u, 17u};

    EXPECT_FALSE(matrix_vector != matrix_vector);
    EXPECT_NE(matrix_vector, matrix_shorter_cols);
    EXPECT_NE(matrix_vector, matrix_shorter_rows);
    EXPECT_NE(matrix_vector, matrix_unequal);
}

TEST_F(trace_matrix, vector)
{
    alignment_trace_matrix matrix{traces, 9u, 17u};

    trace_matrix_test(std::move(matrix));
}

TEST_F(trace_matrix, score_matrix)
{
    alignment_score_matrix score_matrix{scores, 9u, 17u};
    alignment_trace_matrix matrix{database, query, no_config{}, std::move(score_matrix)};

    trace_matrix_test(std::move(matrix));
}

TEST_F(trace_matrix, equal)
{
    // last entry of second row
    std::vector<trace_directions> traces_unequal{traces};
    traces_unequal[2 * 16] = trace_directions::up;

    alignment_trace_matrix matrix_shorter_cols{traces_shorter_cols, 9u, 7u};
    alignment_trace_matrix matrix_shorter_rows{traces_shorter_rows, 4u, 17u};
    alignment_trace_matrix matrix_unequal{std::move(traces_unequal), 9u, 17u};

    alignment_trace_matrix matrix_vector{traces, 9u, 17u};
    alignment_trace_matrix matrix_score_matrix{database, query, no_config{}, alignment_score_matrix{scores, 9u, 17u}};

    EXPECT_EQ(matrix_vector, matrix_vector);
    EXPECT_EQ(matrix_vector, matrix_score_matrix);
    EXPECT_FALSE(matrix_vector == matrix_shorter_cols);
    EXPECT_FALSE(matrix_vector == matrix_shorter_rows);
    EXPECT_FALSE(matrix_vector == matrix_unequal);

    EXPECT_EQ(matrix_score_matrix, matrix_score_matrix);
    EXPECT_EQ(matrix_score_matrix, matrix_vector);
    EXPECT_FALSE(matrix_score_matrix == matrix_shorter_cols);
    EXPECT_FALSE(matrix_score_matrix == matrix_shorter_rows);
    EXPECT_FALSE(matrix_score_matrix == matrix_unequal);
}

TEST_F(trace_matrix, not_equal)
{
    // last entry of second row
    std::vector<trace_directions> traces_unequal{traces};
    traces_unequal[2 * 16] = trace_directions::up;

    alignment_trace_matrix matrix_shorter_cols{traces_shorter_cols, 9u, 7u};
    alignment_trace_matrix matrix_shorter_rows{traces_shorter_rows, 4u, 17u};
    alignment_trace_matrix matrix_unequal{std::move(traces_unequal), 9u, 17u};

    alignment_trace_matrix matrix_vector{traces, 9u, 17u};
    alignment_trace_matrix matrix_score_matrix{database, query, no_config{}, alignment_score_matrix{scores, 9u, 17u}};

    EXPECT_FALSE(matrix_vector != matrix_vector);
    EXPECT_FALSE(matrix_vector != matrix_score_matrix);
    EXPECT_NE(matrix_vector, matrix_shorter_cols);
    EXPECT_NE(matrix_vector, matrix_shorter_rows);
    EXPECT_NE(matrix_vector, matrix_unequal);

    EXPECT_FALSE(matrix_score_matrix != matrix_score_matrix);
    EXPECT_FALSE(matrix_score_matrix != matrix_vector);
    EXPECT_NE(matrix_score_matrix, matrix_shorter_cols);
    EXPECT_NE(matrix_score_matrix, matrix_shorter_rows);
    EXPECT_NE(matrix_score_matrix, matrix_unequal);
}
