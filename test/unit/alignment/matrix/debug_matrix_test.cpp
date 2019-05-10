// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;
using namespace seqan3::detail;

struct debug_matrix_test : public ::testing::Test
{
    std::vector<dna4> sequence1 = "AACACGTTAACCGGTT"_dna4;
    std::vector<dna4> sequence2 = "ACGTACGT"_dna4;
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

    template <typename score_matrix_t>
    void score_matrix_test(score_matrix_t matrix)
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

    template <typename trace_matrix_t>
    void trace_matrix_test(trace_matrix_t matrix)
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

using score_matrix = debug_matrix_test;
using trace_matrix = debug_matrix_test;

template <typename>
struct debug_matrix_traits;

template <Matrix matrix_t, typename sequence1_t, typename sequence2_t>
struct debug_matrix_traits<debug_matrix<matrix_t, sequence1_t, sequence2_t>>
{
    using matrix_type = matrix_t;
    using sequence1_type = sequence1_t;
    using sequence2_type = sequence2_t;
};

TEST_F(debug_matrix_test, matrix_concept)
{
    EXPECT_TRUE((Matrix<row_wise_matrix<int>>));
    EXPECT_TRUE((Matrix<row_wise_matrix<int> &>));
    EXPECT_TRUE((Matrix<row_wise_matrix<int> const>));
    EXPECT_TRUE((Matrix<row_wise_matrix<int> const &>));
    EXPECT_TRUE((Matrix<debug_matrix<row_wise_matrix<int>>>));
    EXPECT_TRUE((Matrix<debug_matrix<row_wise_matrix<int> &>>));
    EXPECT_TRUE((Matrix<debug_matrix<row_wise_matrix<int> const>>));
    EXPECT_TRUE((Matrix<debug_matrix<row_wise_matrix<int> const &>>));
}

TEST_F(debug_matrix_test, construct_with_references)
{
    {
        using debug_matrix_type = decltype(debug_matrix{scores, 9u, 17u, sequence1, sequence2});
        using matrix_type = typename debug_matrix_traits<debug_matrix_type>::matrix_type;
        using sequence1_type = typename debug_matrix_traits<debug_matrix_type>::sequence1_type;
        using sequence2_type = typename debug_matrix_traits<debug_matrix_type>::sequence2_type;

        // the first one has no reference, because {scores, 9u, 17u} will create a row_wise_matrix internally.
        EXPECT_TRUE((std::Same<matrix_type, row_wise_matrix<int>>));
        EXPECT_TRUE((std::Same<sequence1_type, std::vector<dna4> &>));
        EXPECT_TRUE((std::Same<sequence2_type, std::vector<dna4> &>));
    }

    {
        row_wise_matrix row_wise{scores, 9u, 17u};
        using debug_matrix_type = decltype(debug_matrix{row_wise, sequence1, sequence2});
        using matrix_type = typename debug_matrix_traits<debug_matrix_type>::matrix_type;
        using sequence1_type = typename debug_matrix_traits<debug_matrix_type>::sequence1_type;
        using sequence2_type = typename debug_matrix_traits<debug_matrix_type>::sequence2_type;

        EXPECT_TRUE((std::Same<matrix_type, row_wise_matrix<int> &>));
        EXPECT_TRUE((std::Same<sequence1_type, std::vector<dna4> &>));
        EXPECT_TRUE((std::Same<sequence2_type, std::vector<dna4> &>));
    }
}

TEST_F(debug_matrix_test, construct_with_move)
{
    {
        using debug_matrix_type = decltype(debug_matrix{std::move(scores), 9u, 17u,
                                                        std::move(sequence1), std::move(sequence2)});
        using matrix_type = typename debug_matrix_traits<debug_matrix_type>::matrix_type;
        using sequence1_type = typename debug_matrix_traits<debug_matrix_type>::sequence1_type;
        using sequence2_type = typename debug_matrix_traits<debug_matrix_type>::sequence2_type;

        EXPECT_TRUE((std::Same<matrix_type, row_wise_matrix<int>>));
        EXPECT_TRUE((std::Same<sequence1_type, std::vector<dna4>>));
        EXPECT_TRUE((std::Same<sequence2_type, std::vector<dna4>>));
    }

    {
        row_wise_matrix row_wise{scores, 9u, 17u};
        using debug_matrix_type = decltype(debug_matrix{std::move(row_wise),
                                                        std::move(sequence1), std::move(sequence2)});
        using matrix_type = typename debug_matrix_traits<debug_matrix_type>::matrix_type;
        using sequence1_type = typename debug_matrix_traits<debug_matrix_type>::sequence1_type;
        using sequence2_type = typename debug_matrix_traits<debug_matrix_type>::sequence2_type;

        EXPECT_TRUE((std::Same<matrix_type, row_wise_matrix<int>>));
        EXPECT_TRUE((std::Same<sequence1_type, std::vector<dna4>>));
        EXPECT_TRUE((std::Same<sequence2_type, std::vector<dna4>>));
    }
}

TEST_F(score_matrix, vector)
{
    debug_matrix matrix{scores, 9u, 17u};

    score_matrix_test(std::move(matrix));
}

TEST_F(score_matrix, sequences_vector)
{
    debug_matrix matrix{scores, 9u, 17u, sequence1, sequence2};

    EXPECT_EQ(matrix.sequence1(), sequence1);
    EXPECT_EQ(matrix.sequence2(), sequence2);

    score_matrix_test(std::move(matrix));
}

TEST_F(score_matrix, other_matrix)
{
    debug_matrix matrix{scores, 9u, 17u};
    debug_matrix<decltype(matrix)> matrix2{matrix};

    score_matrix_test(std::move(matrix2));
}

TEST_F(score_matrix, sequences_other_matrix)
{
    debug_matrix matrix{scores, 9u, 17u};
    debug_matrix matrix2{matrix, sequence1, sequence2};

    EXPECT_EQ(matrix2.sequence1(), sequence1);
    EXPECT_EQ(matrix2.sequence2(), sequence2);

    score_matrix_test(std::move(matrix));
}

TEST_F(score_matrix, equal)
{
    // last entry of second row
    std::vector<int> scores_unequal{scores};
    scores_unequal[2 * 16] = -16;

    debug_matrix matrix_shorter_cols{scores_shorter_cols, 9u, 7u};
    debug_matrix matrix_shorter_rows{scores_shorter_rows, 4u, 17u};
    debug_matrix matrix_unequal{std::move(scores_unequal), 9u, 17u};

    debug_matrix matrix_vector{scores, 9u, 17u};

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

    debug_matrix matrix_shorter_cols{scores_shorter_cols, 9u, 7u};
    debug_matrix matrix_shorter_rows{scores_shorter_rows, 4u, 17u};
    debug_matrix matrix_unequal{std::move(scores_unequal), 9u, 17u};

    debug_matrix matrix_vector{scores, 9u, 17u};

    EXPECT_FALSE(matrix_vector != matrix_vector);
    EXPECT_NE(matrix_vector, matrix_shorter_cols);
    EXPECT_NE(matrix_vector, matrix_shorter_rows);
    EXPECT_NE(matrix_vector, matrix_unequal);
}

TEST_F(trace_matrix, vector)
{
    debug_matrix matrix{traces, 9u, 17u};

    trace_matrix_test(std::move(matrix));
}

TEST_F(trace_matrix, sequences_vector)
{
    debug_matrix matrix{scores, 9u, 17u, sequence1, sequence2};

    EXPECT_EQ(matrix.sequence1(), sequence1);
    EXPECT_EQ(matrix.sequence2(), sequence2);

    score_matrix_test(std::move(matrix));
}

TEST_F(trace_matrix, other_matrix)
{
    debug_matrix matrix{traces, 9u, 17u};
    debug_matrix<decltype(matrix)> matrix2{traces, 9u, 17u};

    trace_matrix_test(std::move(matrix2));
}

TEST_F(trace_matrix, sequences_other_matrix)
{
    debug_matrix matrix{scores, 9u, 17u};
    debug_matrix matrix2{matrix, sequence1, sequence2};

    EXPECT_EQ(matrix2.sequence1(), sequence1);
    EXPECT_EQ(matrix2.sequence2(), sequence2);

    score_matrix_test(std::move(matrix));
}

TEST_F(trace_matrix, equal)
{
    // last entry of second row
    std::vector<trace_directions> traces_unequal{traces};
    traces_unequal[2 * 16] = trace_directions::up;

    debug_matrix matrix_shorter_cols{traces_shorter_cols, 9u, 7u};
    debug_matrix matrix_shorter_rows{traces_shorter_rows, 4u, 17u};
    debug_matrix matrix_unequal{std::move(traces_unequal), 9u, 17u};

    debug_matrix matrix_vector{traces, 9u, 17u};

    EXPECT_EQ(matrix_vector, matrix_vector);
    EXPECT_FALSE(matrix_vector == matrix_shorter_cols);
    EXPECT_FALSE(matrix_vector == matrix_shorter_rows);
    EXPECT_FALSE(matrix_vector == matrix_unequal);
}

TEST_F(trace_matrix, not_equal)
{
    // last entry of second row
    std::vector<trace_directions> traces_unequal{traces};
    traces_unequal[2 * 16] = trace_directions::up;

    debug_matrix matrix_shorter_cols{traces_shorter_cols, 9u, 7u};
    debug_matrix matrix_shorter_rows{traces_shorter_rows, 4u, 17u};
    debug_matrix matrix_unequal{std::move(traces_unequal), 9u, 17u};

    debug_matrix matrix_vector{traces, 9u, 17u};

    EXPECT_FALSE(matrix_vector != matrix_vector);
    EXPECT_NE(matrix_vector, matrix_shorter_cols);
    EXPECT_NE(matrix_vector, matrix_shorter_rows);
    EXPECT_NE(matrix_vector, matrix_unequal);
}
