// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>

#include "../../../range/iterator_test_template.hpp"

using coordinate_t = seqan3::detail::matrix_coordinate;
using matrix_t = seqan3::detail::coordinate_matrix;
using matrix_iterator_t = std::ranges::iterator_t<matrix_t>;

template <>
struct iterator_fixture<matrix_iterator_t> : public ::testing::Test
{
    using coordinate_column_t = std::vector<coordinate_t>;
    using row_index_t = seqan3::detail::row_index_type<size_t>;
    using column_index_t = seqan3::detail::column_index_type<size_t>;

    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    // Single column with 5 entries as the seq2 has size 4 (need one more for the initialisation row).
    coordinate_column_t column1 = coordinate_column_t{{row_index_t{0u}, column_index_t{0u}},
                                                      {row_index_t{1u}, column_index_t{0u}},
                                                      {row_index_t{2u}, column_index_t{0u}}};
    coordinate_column_t column2 = coordinate_column_t{{row_index_t{0u}, column_index_t{1u}},
                                                      {row_index_t{1u}, column_index_t{1u}},
                                                      {row_index_t{2u}, column_index_t{1u}}};
    std::vector<coordinate_column_t> expected_range{column1, column2};
    matrix_t test_range;

    void SetUp()
    {
        test_range.resize(column_index_t{2u}, row_index_t{3u});
    }

    template <typename actual_column_t, typename expected_column_t>
    static void expect_eq(actual_column_t && actual_column, expected_column_t && expected_column)
    {
        auto actual_it = actual_column.begin();
        auto expected_it = expected_column.begin();
        for (; actual_it != actual_column.end(); ++actual_it, ++expected_it)
        {
            EXPECT_EQ((*actual_it).row, (*expected_it).row);
            EXPECT_EQ((*actual_it).col, (*expected_it).col);
        }
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(matrix_coordinate_test, iterator_fixture, matrix_iterator_t, );

TEST(matrix_coordinate_test, column_concept)
{
    EXPECT_TRUE(std::ranges::forward_range<std::iter_reference_t<matrix_iterator_t>>);
}
