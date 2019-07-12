// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/edit_distance_score_matrix_full.hpp>
#include <seqan3/alignment/matrix/matrix_concept.hpp>

using namespace seqan3;

using score_type = int;
using word_type = uint8_t;

template <bool is_semi_global, bool use_max_errors>
class matrix_type
    : public seqan3::detail::edit_distance_score_matrix_full<word_type, score_type, is_semi_global, use_max_errors>
{
public:
    using base_t = seqan3::detail::edit_distance_score_matrix_full<word_type, score_type,
                                                                   is_semi_global, use_max_errors>;

    matrix_type(size_t const rows_size) : base_t{rows_size}
    {}

    using base_t::max_rows;
    using base_t::add_column;
    using base_t::reserve;
};

static constexpr score_type INF = detail::matrix_inf<score_type>;

std::vector<std::vector<int>> as_row_wise_vector(auto matrix)
{
    std::vector<std::vector<int>> result{};
    for (unsigned row = 0; row < matrix.rows(); ++row)
    {
        result.push_back({});
        for (unsigned col = 0; col < matrix.cols(); ++col)
        {
            std::optional<int> entry = matrix.at(row, col);
            result.back().push_back(entry.value_or(INF));
        }
    }
    return result;
}

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
    std::vector<std::vector<int>> expect
    {
        {-0, -1, -2, -3, -4, -5, -6, -7, -8, -9},
        {-1, -0, -1, -2, -3, -4, -5, -6, -7, -8},
        {-2, -1, -1, -1, -2, -3, -4, -5, -6, -7},
        {-3, -2, -2, -2, -2, -2, -3, -4, -5, -6},
        {-4, -3, -3, -3, -3, -3, -3, -3, -4, -5},
        {-5, -4, -3, -4, -4, -4, -4, -4, -4, -4},
        {-6, -5, -4, -3, -4, -5, -5, -5, -5, -5},
        {-7, -6, -5, -4, -4, -4, -5, -6, -6, -6},
        {-8, -7, -6, -5, -5, -5, -5, -5, -6, -7}
    };

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
    std::vector<std::vector<int>> expect
    {
        { -0,  -1,  -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9},
        { -1,  -0,  -1,  -2,  -3,  -4,  -5,  -6,  -7,  -8},
        { -2,  -1,  -1,  -2,  -3,  -3,  -4,  -5,  -6,  -7},
        { -3,  -2,  -1,  -2,  -3,  -4,  -3,  -4,  -5,  -6},
        { -4,  -3,  -2,  -2,  -3,  -4,  -4,  -4,  -5,  -6},
        { -5,  -4,  -3,  -2,  -3,  -4,  -5,  -4,  -5,  -6},
        { -6,  -5,  -4,  -3,  -3,  -4,  -5,  -5,  -5,  -6},
        { -7,  -6,  -5,  -4,  -3,  -4,  -5,  -6,  -5,  -6},
        { -8,  -7,  -6,  -5,  -4,  -4,  -5,  -6,  -6,  -6},
        { -9,  -8,  -7,  -6,  -5,  -4,  -5,  -6,  -7,  -6},
        {-10,  -9,  -8,  -7,  -6,  -5,  -5,  -6,  -7,  -7},
        {-11, -10,  -9,  -8,  -7,  -6,  -5,  -6,  -7,  -8},
        {-12, -11, -10,  -9,  -8,  -7,  -6,  -6,  -7,  -8},
        {-13, -12, -11, -10,  -9,  -8,  -7,  -6,  -7,  -8},
        {-14, -13, -12, -11, -10,  -9,  -8,  -7,  -7,  -8},
        {-15, -14, -13, -12, -11, -10,  -9,  -8,  -7,  -8},
        {-16, -15, -14, -13, -12, -11, -10,  -9,  -8,  -8},
        {-17, -16, -15, -14, -13, -12, -11, -10,  -9,  -8}
    };

    EXPECT_EQ(result, expect);
}

TEST(semi_global, empty)
{
    matrix_type<true, false> matrix{1u};

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global, epsilon)
{
    matrix_type<true, false> matrix{1u};

    matrix.add_column({}, {});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global, epsilon_row)
{
    matrix_type<true, false> matrix{1u};

    matrix.add_column({}, {});
    matrix.add_column({}, {});
    matrix.add_column({}, {});
    matrix.add_column({}, {});
    matrix.add_column({}, {});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect{{-0, -0, -0, -0, -0}};

    EXPECT_EQ(result, expect);
}

TEST(semi_global, single_word)
{
    matrix_type<true, false> matrix{9u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u}, {0b0000'0000u});
    matrix.add_column({0b1111'1110u}, {0b0000'0000u});
    matrix.add_column({0b1110'1110u}, {0b0000'0000u});
    matrix.add_column({0b1101'1101u}, {0b0000'0010u});
    matrix.add_column({0b1101'1001u}, {0b0000'0000u});
    matrix.add_column({0b1011'1011u}, {0b0100'0100u});
    matrix.add_column({0b0011'0011u}, {0b0000'0000u});
    matrix.add_column({0b0111'0111u}, {0b1000'1000u});
    matrix.add_column({0b0110'0111u}, {0b0000'0000u});
    matrix.add_column({0b1110'1110u}, {0b0000'0000u});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        {-0, -0, -0, -0, -0, -0, -0, -0, -0, -0},
        {-1, -0, -0, -1, -1, -1, -1, -1, -1, -0},
        {-2, -1, -1, -0, -1, -2, -2, -2, -2, -1},
        {-3, -2, -2, -1, -1, -1, -2, -3, -3, -2},
        {-4, -3, -3, -2, -2, -2, -2, -2, -3, -3},
        {-5, -4, -3, -3, -3, -3, -3, -3, -3, -3},
        {-6, -5, -4, -3, -3, -4, -4, -4, -4, -4},
        {-7, -6, -5, -4, -4, -3, -4, -5, -5, -5},
        {-8, -7, -6, -5, -5, -4, -4, -4, -5, -6}
    };

    EXPECT_EQ(result, expect);
}

TEST(semi_global, multiple_words)
{
    matrix_type<true, false> matrix{18u};
    matrix.reserve(10u);

    matrix.add_column({0b1111'1111u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1111'1110u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1111'1001u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1110'0011u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1000'0111u, 0b1111'1111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0001'1110u, 0b1111'1110u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0111'1101u, 0b1111'1000u, 0b1u}, {0b0000'0010u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1111'0001u, 0b1110'0001u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b1100'0011u, 0b1000'0111u, 0b1u}, {0b0000'0000u, 0b0000'0000u, 0b0u});
    matrix.add_column({0b0100'1110u, 0b0001'1111u, 0b0u}, {0b0001'0000u, 0b0000'0000u, 0b0u});

    // row-wise matrix
    std::vector<std::vector<int>> result = as_row_wise_vector(matrix);
    std::vector<std::vector<int>> expect
    {
        { -0,  -0,  -0,  -0,  -0,  -0,  -0,  -0,  -0,  -0},
        { -1,  -0,  -1,  -1,  -1,  -0,  -1,  -1,  -1,  -0},
        { -2,  -1,  -1,  -2,  -2,  -1,  -0,  -1,  -2,  -1},
        { -3,  -2,  -1,  -2,  -3,  -2,  -1,  -1,  -2,  -2},
        { -4,  -3,  -2,  -2,  -3,  -3,  -2,  -1,  -2,  -3},
        { -5,  -4,  -3,  -2,  -3,  -4,  -3,  -2,  -2,  -2},
        { -6,  -5,  -4,  -3,  -3,  -4,  -4,  -3,  -2,  -2},
        { -7,  -6,  -5,  -4,  -3,  -4,  -5,  -4,  -3,  -3},
        { -8,  -7,  -6,  -5,  -4,  -4,  -5,  -5,  -4,  -3},
        { -9,  -8,  -7,  -6,  -5,  -4,  -5,  -6,  -5,  -4},
        {-10,  -9,  -8,  -7,  -6,  -5,  -5,  -6,  -6,  -5},
        {-11, -10,  -9,  -8,  -7,  -6,  -5,  -6,  -7,  -6},
        {-12, -11, -10,  -9,  -8,  -7,  -6,  -6,  -7,  -7},
        {-13, -12, -11, -10,  -9,  -8,  -7,  -6,  -7,  -8},
        {-14, -13, -12, -11, -10,  -9,  -8,  -7,  -7,  -8},
        {-15, -14, -13, -12, -11, -10,  -9,  -8,  -7,  -8},
        {-16, -15, -14, -13, -12, -11, -10,  -9,  -8,  -8},
        {-17, -16, -15, -14, -13, -12, -11, -10,  -9,  -8}
    };

    EXPECT_EQ(result, expect);
}

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
    std::vector<std::vector<int>> expect
    {
        { -0, -0, -0, -0, -0, -0, -0, -0, -0, -0},
        { -1, -0, -0, -1, -1, -1, -1, -1, -1, -0},
        { -2, -1, -1, -0, -1, -2, -2, -2, -2, -1},
        { -3, -2, -2, -1, -1, -1, -2, -3, -3, -2},
        { -4, -3, -3, -2, -2, -2, -2, -2, -3, -3},
        { -5, -4, -3, -3, -3, -3, -3, -3, -3, -3},
        {INF, -5, -4, -3, -3, -4, -4, -4, -4, -4},
        {INF,INF, -5, -4, -4, -3, -4, -5, -5, -5},
        {INF,INF,INF, -5, -5, -4, -4, -4, -5,INF}
    };

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
    std::vector<std::vector<int>> expect
    {
        { -0, -0, -0, -0, -0, -0, -0, -0, -0, -0},
        { -1, -0, -1, -1, -1, -0, -1, -1, -1, -0},
        { -2, -1, -1, -2, -2, -1, -0, -1, -2, -1},
        { -3, -2, -1, -2, -3, -2, -1, -1, -2, -2},
        { -4, -3, -2, -2, -3, -3, -2, -1, -2, -3},
        { -5, -4, -3, -2, -3, -4, -3, -2, -2, -2},
        { -6, -5, -4, -3, -3, -4, -4, -3, -2, -2},
        { -7, -6, -5, -4, -3, -4, -5, -4, -3, -3},
        { -8, -7, -6, -5, -4, -4, -5, -5, -4, -3},
        {INF, -8, -7, -6, -5, -4, -5, -6, -5, -4},
        {INF,INF, -8, -7, -6, -5, -5, -6, -6, -5},
        {INF,INF,INF, -8, -7, -6, -5, -6, -7, -6},
        {INF,INF,INF,INF, -8, -7, -6, -6, -7, -7},
        {INF,INF,INF,INF,INF, -8, -7, -6, -7, -8},
        {INF,INF,INF,INF,INF,INF, -8, -7, -7, -8},
        {INF,INF,INF,INF,INF,INF,INF, -8, -7, -8},
        {INF,INF,INF,INF,INF,INF,INF,INF, -8, -8},
        {INF,INF,INF,INF,INF,INF,INF,INF,INF, -8}
    };

    EXPECT_EQ(result, expect);
}
