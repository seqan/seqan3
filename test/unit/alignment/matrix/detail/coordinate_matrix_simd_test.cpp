// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/test/simd_utility.hpp>
#include <seqan3/utility/simd/simd.hpp>

#include "../../../range/iterator_test_template.hpp"

using simd_t = seqan3::simd::simd_type_t<int32_t>;
using coordinate_t = seqan3::detail::simd_matrix_coordinate<simd_t>;
using matrix_t = seqan3::detail::coordinate_matrix<simd_t>;
using matrix_iterator_t = std::ranges::iterator_t<matrix_t>;

template <>
struct iterator_fixture<matrix_iterator_t> : public ::testing::Test
{
    using index_column_t = std::vector<coordinate_t, seqan3::aligned_allocator<coordinate_t, alignof(coordinate_t)>>;
    using row_index_t = seqan3::detail::row_index_type<int32_t>;
    using column_index_t = seqan3::detail::column_index_type<int32_t>;

    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    // Single column with 5 entries as the seq2 has size 4 (need one more for the initialisation row).
    index_column_t column1 = index_column_t{{row_index_t{0}, column_index_t{0}},
                                            {row_index_t{1}, column_index_t{0}},
                                            {row_index_t{2}, column_index_t{0}}};
    index_column_t column2 = index_column_t{{row_index_t{0}, column_index_t{1}},
                                            {row_index_t{1}, column_index_t{1}},
                                            {row_index_t{2}, column_index_t{1}}};
    std::vector<index_column_t> expected_range{column1, column2};
    matrix_t test_range;

    void SetUp()
    {
        test_range.resize(column_index_t{2}, row_index_t{3});
    }

    template <typename actual_column_t, typename expected_column_t>
    static void expect_eq(actual_column_t && actual_column, expected_column_t && expected_column)
    {
        auto actual_it = actual_column.begin();
        auto expected_it = expected_column.begin();
        for (; actual_it != actual_column.end(); ++actual_it, ++expected_it)
        {
            SIMD_EQ((*actual_it).row, (*expected_it).row);
            SIMD_EQ((*actual_it).col, (*expected_it).col);
        }
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(matrix_coordinate_test, iterator_fixture, matrix_iterator_t, );

TEST(matrix_coordinate_test, column_concept)
{
    EXPECT_TRUE(std::ranges::forward_range<std::iter_reference_t<matrix_iterator_t>>);
}
