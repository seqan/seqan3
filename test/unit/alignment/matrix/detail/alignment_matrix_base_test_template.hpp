// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/std/ranges>

template <typename t>
struct alignment_matrix_base_test : public ::testing::Test
{
    using matrix_t = std::tuple_element_t<0, t>;
    static constexpr bool is_banded = std::tuple_element_t<1, t>::value;

    alignment_matrix_base_test()
    {
        if constexpr (is_banded)
            test_range = matrix_t{first, second, seqan3::static_band{seqan3::lower_bound{-2}, seqan3::upper_bound{2}}};
        else
            test_range = matrix_t{first, second};
    }

    std::string first{"abba"};
    std::string second{"baba"};
    matrix_t test_range{};
};

TYPED_TEST_SUITE_P(alignment_matrix_base_test);

TYPED_TEST_P(alignment_matrix_base_test, range_concepts)
{
    using outer_it = std::ranges::iterator_t<typename TestFixture::matrix_t>;
    using column_t = std::iter_value_t<outer_it>;
    using inner_it = std::ranges::iterator_t<column_t>;

    EXPECT_TRUE(std::ranges::input_range<typename TestFixture::matrix_t>);
    EXPECT_TRUE(std::input_iterator<outer_it>);
    EXPECT_TRUE(std::input_iterator<inner_it>);
    EXPECT_TRUE(std::ranges::input_range<column_t>);
    EXPECT_TRUE(std::ranges::forward_range<column_t>);
    EXPECT_TRUE(std::ranges::sized_range<column_t>);
    EXPECT_TRUE(std::ranges::view<column_t>);
}

TYPED_TEST_P(alignment_matrix_base_test, begin_end)
{
    auto it = this->test_range.begin();
    auto col = *it;
    auto col_it = col.begin();

    EXPECT_NE(col_it, col.end());
    EXPECT_NE(it, this->test_range.end());
}

TYPED_TEST_P(alignment_matrix_base_test, basic_construction)
{
    using matrix_type = typename TestFixture::matrix_t;

    EXPECT_TRUE(std::is_default_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_copy_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_move_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_copy_assignable_v<matrix_type>);
    EXPECT_TRUE(std::is_move_assignable_v<matrix_type>);
    EXPECT_TRUE(std::is_destructible_v<matrix_type>);
}

template <typename matrix_t>
void test_matrix_iteration(matrix_t test_range, int32_t const expected_col_count, int32_t const expected_row_count)
{
    int32_t col_count = 0;
    int32_t row_count = 0;
    for (auto col : test_range)
    {
        for ([[maybe_unused]] auto cell : col)
            ++row_count;

        ++col_count;
    }
    EXPECT_EQ(col_count, expected_col_count);
    EXPECT_EQ(row_count, expected_row_count);
}

TYPED_TEST_P(alignment_matrix_base_test, empty_row)
{
    using matrix_type = typename TestFixture::matrix_t;

    if constexpr (TestFixture::is_banded)
    {
        seqan3::static_band band{seqan3::lower_bound{-2}, seqan3::upper_bound{4}};
        test_matrix_iteration(matrix_type{this->first, std::string{""}, band}, 5, 5);
    }
    else
    {
        test_matrix_iteration(matrix_type{this->first, std::string{""}}, 5, 5);
    }
}

TYPED_TEST_P(alignment_matrix_base_test, empty_col)
{
    using matrix_type = typename TestFixture::matrix_t;

    if constexpr (TestFixture::is_banded)
    {
        seqan3::static_band band{seqan3::lower_bound{-2}, seqan3::upper_bound{2}};
        test_matrix_iteration(matrix_type{std::string{""}, this->second, band}, 1, 3);
    }
    else
    {
        test_matrix_iteration(matrix_type{std::string{""}, this->second}, 1, 5);
    }
}

TYPED_TEST_P(alignment_matrix_base_test, empty_col_row)
{
    using matrix_type = typename TestFixture::matrix_t;

    if constexpr (TestFixture::is_banded)
    {
        seqan3::static_band band{seqan3::lower_bound{0}, seqan3::upper_bound{2}};
        test_matrix_iteration(matrix_type{std::string{""}, std::string{""}, band}, 1, 1);
    }
    else
    {
        test_matrix_iteration(matrix_type{std::string{""}, std::string{""}}, 1, 1);
    }
}

REGISTER_TYPED_TEST_SUITE_P(alignment_matrix_base_test, range_concepts, begin_end, basic_construction, empty_row,
                            empty_col, empty_col_row);
