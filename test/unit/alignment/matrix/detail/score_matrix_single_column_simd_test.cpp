// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/test/simd_utility.hpp>
#include <seqan3/utility/simd/simd.hpp>

#include "../../../range/iterator_test_template.hpp"

using score_t = seqan3::simd::simd_type_t<int32_t>;
using matrix_t = seqan3::detail::score_matrix_single_column<score_t>;
using matrix_iterator_t = std::ranges::iterator_t<matrix_t>;

template <>
struct iterator_fixture<matrix_iterator_t> : public ::testing::Test
{
    using cell_t = std::tuple<score_t, score_t, score_t>;
    using score_column_t = std::vector<cell_t, seqan3::aligned_allocator<cell_t, alignof(cell_t)>>;

    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    // Single column with 5 entries as the seq2 has size 4 (need one more for the initialisation row).
    score_column_t column = score_column_t{5, std::tuple{score_t{0}, score_t{0}, score_t{0}}};
    std::vector<score_column_t> expected_range{column, column, column, column};
    matrix_t test_range;

    void SetUp()
    {
        test_range.resize(seqan3::detail::column_index_type{4}, seqan3::detail::row_index_type{5});
    }

    template <typename actual_column_t, typename expected_column_t>
    static void expect_eq(actual_column_t && actual_column, expected_column_t && expected_column)
    {
        auto actual_it = actual_column.begin();
        auto expected_it = expected_column.begin();
        for (; actual_it != actual_column.end(); ++actual_it, ++expected_it)
        {
            auto actual_cell = *actual_it;
            auto expected_cell = *expected_it;

            SIMD_EQ(actual_cell.best_score(), std::get<0>(expected_cell));
            SIMD_EQ(actual_cell.horizontal_score(), std::get<1>(expected_cell));
            SIMD_EQ(actual_cell.vertical_score(), std::get<2>(expected_cell));
        }
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(score_matrix_single_column_simd_test, iterator_fixture, matrix_iterator_t, );
