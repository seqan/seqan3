// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/alignment/matrix/detail/combined_score_and_trace_matrix.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/matrix/detail/trace_matrix_full.hpp>

#include "../../../range/iterator_test_template.hpp"

using score_t = int32_t;
using score_matrix_t = seqan3::detail::score_matrix_single_column<score_t>;
using trace_t = seqan3::detail::trace_directions;
using trace_matrix_t = seqan3::detail::trace_matrix_full<trace_t>;

using matrix_t = seqan3::detail::combined_score_and_trace_matrix<score_matrix_t, trace_matrix_t>;
using matrix_iterator_t = std::ranges::iterator_t<matrix_t>;

template <>
struct iterator_fixture<matrix_iterator_t> : public ::testing::Test
{
    using alignment_column_t =
        std::vector<std::pair<std::tuple<score_t, score_t, score_t>, std::tuple<trace_t, trace_t, trace_t>>>;

    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    static constexpr seqan3::detail::trace_directions none = seqan3::detail::trace_directions::none;

    // Single column with 5 entries as the seq2 has size 4 (need one more for the initialisation row).
    alignment_column_t column = alignment_column_t{5, std::pair{std::tuple{0, 0, 0}, std::tuple{none, none, none}}};
    std::vector<alignment_column_t> expected_range{column, column, column, column};
    matrix_t test_range;

    void SetUp()
    {
        std::string seq1 = "abc";
        std::string seq2 = "abcd";

        test_range.resize(seqan3::detail::column_index_type<size_t>{4}, seqan3::detail::row_index_type<size_t>{5});
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

            EXPECT_EQ(actual_cell.best_score(), std::get<0>(std::get<0>(expected_cell)));
            EXPECT_EQ(actual_cell.horizontal_score(), std::get<1>(std::get<0>(expected_cell)));
            EXPECT_EQ(actual_cell.vertical_score(), std::get<2>(std::get<0>(expected_cell)));
            EXPECT_EQ(actual_cell.best_trace(), std::get<0>(std::get<1>(expected_cell)));
            EXPECT_EQ(actual_cell.horizontal_trace(), std::get<1>(std::get<1>(expected_cell)));
            EXPECT_EQ(actual_cell.vertical_trace(), std::get<2>(std::get<1>(expected_cell)));
        }
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(trace_matrix_single_column_test, iterator_fixture, matrix_iterator_t, );
