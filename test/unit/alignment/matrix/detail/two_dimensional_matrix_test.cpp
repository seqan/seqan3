// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <vector>

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_concept.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/matrix_concept.hpp>

#include "../../../range/iterator_test_template.hpp"

using namespace seqan3;
using namespace seqan3::detail;

using row_wise = std::integral_constant<matrix_major_order, matrix_major_order::row>;
using col_wise = std::integral_constant<matrix_major_order, matrix_major_order::column>;

template <typename constant_t>
struct two_dimensional_matrix_test : public ::testing::Test
{
    using matrix_type = seqan3::detail::two_dimensional_matrix<int, std::allocator<int>, constant_t::value>;

    std::vector<int> expected_range{0,  1,  2,  3,
                                    4,  5,  6,  7,
                                    8,  9, 10, 11};
};

using testing_types = ::testing::Types<row_wise, col_wise>;
TYPED_TEST_CASE(two_dimensional_matrix_test, testing_types);

TYPED_TEST(two_dimensional_matrix_test, concepts)
{
    using matrix_type = typename TestFixture::matrix_type;

    EXPECT_TRUE(matrix<matrix_type>);
    EXPECT_TRUE(std::ranges::input_range<matrix_type>);
    EXPECT_TRUE(std::ranges::forward_range<matrix_type>);
    EXPECT_TRUE(std::ranges::bidirectional_range<matrix_type>);
    EXPECT_TRUE(std::ranges::random_access_range<matrix_type>);
}

TYPED_TEST(two_dimensional_matrix_test, construction)
{
    using matrix_type = typename TestFixture::matrix_type;

    EXPECT_TRUE(std::is_default_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_copy_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_move_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_copy_assignable_v<matrix_type>);
    EXPECT_TRUE(std::is_move_assignable_v<matrix_type>);
    EXPECT_TRUE(std::is_destructible_v<matrix_type>);

    EXPECT_TRUE((std::is_constructible_v<matrix_type, number_rows, number_cols>));
}

TYPED_TEST(two_dimensional_matrix_test, cols)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{number_rows{3}, number_cols{4}};
    EXPECT_EQ(matrix.cols(), 4u);
}

TYPED_TEST(two_dimensional_matrix_test, rows)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{number_rows{3}, number_cols{4}};
    EXPECT_EQ(matrix.rows(), 3u);
}

TYPED_TEST(two_dimensional_matrix_test, range)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{number_rows{3}, number_cols{4}};
    auto it = std::copy(this->expected_range.begin(), this->expected_range.end(), matrix.begin());

    EXPECT_EQ(it, matrix.end());
    EXPECT_TRUE(std::equal(this->expected_range.begin(), this->expected_range.end(), matrix.begin()));
}

TYPED_TEST(two_dimensional_matrix_test, subscript)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{number_rows{3}, number_cols{4}};
    auto it = std::copy(this->expected_range.begin(), this->expected_range.end(), matrix.begin());

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ((matrix[{row_index_type{0u}, column_index_type{0u}}]), 0);
        EXPECT_EQ((matrix[{row_index_type{1u}, column_index_type{0u}}]), 4);
        EXPECT_EQ((matrix[{row_index_type{0u}, column_index_type{1u}}]), 1);
        EXPECT_EQ((matrix[{row_index_type{0u}, column_index_type{3u}}]), 3);
        EXPECT_EQ((matrix[{row_index_type{2u}, column_index_type{0u}}]), 8);
        EXPECT_EQ((matrix[{row_index_type{2u}, column_index_type{3u}}]), 11);
    }
    else  // Column order.
    {
        EXPECT_EQ((matrix[{row_index_type{0u}, column_index_type{0u}}]), 0);
        EXPECT_EQ((matrix[{row_index_type{1u}, column_index_type{0u}}]), 1);
        EXPECT_EQ((matrix[{row_index_type{0u}, column_index_type{1u}}]), 3);
        EXPECT_EQ((matrix[{row_index_type{0u}, column_index_type{3u}}]), 9);
        EXPECT_EQ((matrix[{row_index_type{2u}, column_index_type{0u}}]), 2);
        EXPECT_EQ((matrix[{row_index_type{2u}, column_index_type{3u}}]), 11);
    }
}

TYPED_TEST(two_dimensional_matrix_test, at)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{number_rows{3}, number_cols{4}};
    auto it = std::copy(this->expected_range.begin(), this->expected_range.end(), matrix.begin());

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ((matrix.at({row_index_type{0u}, column_index_type{0u}})), 0);
        EXPECT_EQ((matrix.at({row_index_type{1u}, column_index_type{0u}})), 4);
        EXPECT_EQ((matrix.at({row_index_type{0u}, column_index_type{1u}})), 1);
        EXPECT_EQ((matrix.at({row_index_type{0u}, column_index_type{3u}})), 3);
        EXPECT_EQ((matrix.at({row_index_type{2u}, column_index_type{0u}})), 8);
        EXPECT_EQ((matrix.at({row_index_type{2u}, column_index_type{3u}})), 11);
    }
    else
    {
        EXPECT_EQ((matrix.at({row_index_type{0u}, column_index_type{0u}})), 0);
        EXPECT_EQ((matrix.at({row_index_type{1u}, column_index_type{0u}})), 1);
        EXPECT_EQ((matrix.at({row_index_type{0u}, column_index_type{1u}})), 3);
        EXPECT_EQ((matrix.at({row_index_type{0u}, column_index_type{3u}})), 9);
        EXPECT_EQ((matrix.at({row_index_type{2u}, column_index_type{0u}})), 2);
        EXPECT_EQ((matrix.at({row_index_type{2u}, column_index_type{3u}})), 11);

    }

    EXPECT_THROW((matrix.at({row_index_type{3u}, column_index_type{3u}})), std::invalid_argument);
    EXPECT_THROW((matrix.at({row_index_type{2u}, column_index_type{4u}})), std::invalid_argument);
}

TYPED_TEST(two_dimensional_matrix_test, conversion)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{number_rows{3}, number_cols{4}};
    auto it = std::copy(this->expected_range.begin(), this->expected_range.end(), matrix.begin());

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        using converted_t = two_dimensional_matrix<uint32_t, std::allocator<uint32_t>, matrix_major_order::column>;
        converted_t converted = static_cast<converted_t>(matrix);
        EXPECT_EQ(converted.rows(), matrix.rows());
        EXPECT_EQ(converted.cols(), matrix.cols());
        EXPECT_EQ((converted[{row_index_type{0u}, column_index_type{0u}}]), 0u);
        EXPECT_EQ((converted[{row_index_type{1u}, column_index_type{0u}}]), 4u);
        EXPECT_EQ((converted[{row_index_type{0u}, column_index_type{1u}}]), 1u);
        EXPECT_EQ((converted[{row_index_type{0u}, column_index_type{3u}}]), 3u);
        EXPECT_EQ((converted[{row_index_type{2u}, column_index_type{0u}}]), 8u);
        EXPECT_EQ((converted[{row_index_type{2u}, column_index_type{3u}}]), 11u);
    }
    else
    {
        using converted_t = two_dimensional_matrix<uint32_t, std::allocator<uint32_t>, matrix_major_order::row>;
        converted_t converted = static_cast<converted_t>(matrix);
        EXPECT_EQ(converted.rows(), matrix.rows());
        EXPECT_EQ(converted.cols(), matrix.cols());
        EXPECT_EQ((converted[{row_index_type{0u}, column_index_type{0u}}]), 0u);
        EXPECT_EQ((converted[{row_index_type{1u}, column_index_type{0u}}]), 1u);
        EXPECT_EQ((converted[{row_index_type{0u}, column_index_type{1u}}]), 3u);
        EXPECT_EQ((converted[{row_index_type{0u}, column_index_type{3u}}]), 9u);
        EXPECT_EQ((converted[{row_index_type{2u}, column_index_type{0u}}]), 2u);
        EXPECT_EQ((converted[{row_index_type{2u}, column_index_type{3u}}]), 11u);

    }
}

//-----------------------------------------------------------------------------
// Iterator tests
//-----------------------------------------------------------------------------

template <typename policy_t>
struct iterator_fixture<two_dimensional_matrix_test<policy_t>> : two_dimensional_matrix_test<policy_t>
{
    using base_t = two_dimensional_matrix_test<policy_t>;
    using matrix_type = typename base_t::matrix_type;

    // Test random access range.
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    matrix_type test_range{number_rows{3}, number_cols{4}, base_t::expected_range};
};

using iter_testing_types = ::testing::Types<two_dimensional_matrix_test<row_wise>,
                                            two_dimensional_matrix_test<col_wise>>;

INSTANTIATE_TYPED_TEST_CASE_P(two_dimensional_iterator, iterator_fixture, iter_testing_types);

template <typename policy_t>
struct iterator_fixture_two_dimensional : two_dimensional_matrix_test<policy_t>
{
    using base_t = two_dimensional_matrix_test<policy_t>;
    using matrix_type = typename base_t::matrix_type;
    using iterator_type = typename matrix_type::iterator;

    matrix_type test_range{number_rows{3}, number_cols{4}, base_t::expected_range};
};

TYPED_TEST_CASE(iterator_fixture_two_dimensional, testing_types);

TYPED_TEST(iterator_fixture_two_dimensional, two_dimensional_concept)
{
    EXPECT_TRUE(two_dimensional_matrix_iterator<typename TestFixture::iterator_type>);
}

TYPED_TEST(iterator_fixture_two_dimensional, update_by_matrix_offset_add)
{
    auto it = this->test_range.begin();
    auto it_advanced = it += matrix_offset{row_index_type{1}, column_index_type{2}};

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ(*it, 6);
        EXPECT_EQ(*it_advanced, 6);
    }
    else
    {
        EXPECT_EQ(*it, 7);
        EXPECT_EQ(*it_advanced, 7);
    }
}

TYPED_TEST(iterator_fixture_two_dimensional, advance_by_matrix_offset_add)
{
    auto it = this->test_range.begin();
    auto it_advanced = it + matrix_offset{row_index_type{1}, column_index_type{2}};

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ(*it, 0);
        EXPECT_EQ(*it_advanced, 6);
    }
    else
    {
        EXPECT_EQ(*it, 0);
        EXPECT_EQ(*it_advanced, 7);
    }
}

TYPED_TEST(iterator_fixture_two_dimensional, advance_by_matrix_offset_add_friend)
{
    auto it = this->test_range.begin();
    auto it_advanced = matrix_offset{row_index_type{1}, column_index_type{2}} + it;

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ(*it, 0);
        EXPECT_EQ(*it_advanced, 6);
    }
    else
    {
        EXPECT_EQ(*it, 0);
        EXPECT_EQ(*it_advanced, 7);
    }
}

TYPED_TEST(iterator_fixture_two_dimensional, update_by_matrix_offset_subtract)
{
    auto it = this->test_range.begin() + matrix_offset{row_index_type{2}, column_index_type{3}};
    auto it_advanced = it -= matrix_offset{row_index_type{1}, column_index_type{2}};

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ(*it, 5);
        EXPECT_EQ(*it_advanced, 5);
    }
    else
    {
        EXPECT_EQ(*it, 4);
        EXPECT_EQ(*it_advanced, 4);
    }
}

TYPED_TEST(iterator_fixture_two_dimensional, advance_by_matrix_offset_subtract)
{
    auto it = this->test_range.begin() + matrix_offset{row_index_type{2}, column_index_type{3}};
    auto it_advanced = it - matrix_offset{row_index_type{1}, column_index_type{2}};

    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        EXPECT_EQ(*it, 11);
        EXPECT_EQ(*it_advanced, 5);
    }
    else
    {
        EXPECT_EQ(*it, 11);
        EXPECT_EQ(*it_advanced, 4);
    }
}

TYPED_TEST(iterator_fixture_two_dimensional, coordinate)
{
    auto it = this->test_range.begin();
    std::ptrdiff_t pos = 0;
    if constexpr (TypeParam::value == matrix_major_order::row)
    {
        for (; it != this->test_range.end(); ++it, ++pos)
        {
            EXPECT_EQ(it.coordinate().col, pos % this->test_range.cols());
            EXPECT_EQ(it.coordinate().row, pos / this->test_range.cols());
        }
    }
    else
    {
        for (; it != this->test_range.end(); ++it, ++pos)
        {
            EXPECT_EQ(it.coordinate().col, pos / this->test_range.rows());
            EXPECT_EQ(it.coordinate().row, pos % this->test_range.rows());
        }
    }
}
