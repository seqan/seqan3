// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include <seqan3/alignment/matrix/detail/two_dimensional_matrix_iterator_concept.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/test/simd_utility.hpp>

#include "../../../range/iterator_test_template.hpp"

template <typename score_t, seqan3::detail::matrix_major_order order, typename allocator_t = std::allocator<score_t>>
using test_matrix_t = seqan3::detail::two_dimensional_matrix<score_t, allocator_t, order>;

template <typename score_type, seqan3::detail::matrix_major_order order = seqan3::detail::matrix_major_order::row>
std::vector<score_type, std::allocator<score_type>> create_matrix_storage()
{
    // this is a small hack to allow simd types in an initialiser list on gcc 7 and gcc 10;
    using storage_t = std::array<score_type, 12>;

    // note: we represent the same matrix in one case with a row-wise data layout and in the other case with a
    // column-wise data layout. Also note that the score_type can be a builtin integer or a simd vector of a builtin
    // integer.
    storage_t row_wise
    {{
        score_type{0}, score_type{1}, score_type{ 2}, score_type{ 3},
        score_type{4}, score_type{5}, score_type{ 6}, score_type{ 7},
        score_type{8}, score_type{9}, score_type{10}, score_type{11}
    }};

    storage_t column_wise
    {{
        score_type{0}, score_type{4}, score_type{8},
        score_type{1}, score_type{5}, score_type{9},
        score_type{2}, score_type{6}, score_type{10},
        score_type{3}, score_type{7}, score_type{11}
    }};

    // for a simd vector we make sure that some simd values are completely set, i.e. each scalar value in that simd
    // vector has some value.
    if constexpr(seqan3::simd::simd_concept<score_type>)
    {
        row_wise[5] = seqan3::simd::iota<score_type>(5); // 5, 6, 7, ...
        row_wise[8] = seqan3::simd::fill<score_type>(8); // 8, 8, 8, ...

        column_wise[4] = seqan3::simd::iota<score_type>(5);
        column_wise[2] = seqan3::simd::fill<score_type>(8);
    }

    storage_t storage = order == seqan3::detail::matrix_major_order::row ? row_wise : column_wise;
    return {storage.begin(), storage.end()};
}

template <typename score_t>
struct make_unsigned_score_type : std::make_unsigned<score_t> {};

template <seqan3::simd::simd_concept simd_score_t>
struct make_unsigned_score_type<simd_score_t>
{
    using score_type = typename seqan3::simd::simd_traits<simd_score_t>::scalar_type;
    using unsigned_score_type = std::make_unsigned_t<score_type>;
    static constexpr auto length = seqan3::simd::simd_traits<simd_score_t>::length;

    using type = seqan3::simd::simd_type_t<unsigned_score_type, length>;
};

//-----------------------------------------------------------------------------
// Matrix tests
//-----------------------------------------------------------------------------

template <typename matrix_t>
struct two_dimensional_matrix_test;

template <typename score_t, seqan3::detail::matrix_major_order order>
struct two_dimensional_matrix_test<test_matrix_t<score_t, order>> : public ::testing::Test
{
    using matrix_type = test_matrix_t<score_t, order>;
    using score_type = typename matrix_type::value_type;
    static constexpr seqan3::detail::matrix_major_order matrix_order = order;

    std::vector<score_type> expected_matrix_content{create_matrix_storage<score_type,
                                                    seqan3::detail::matrix_major_order::row>()};
    std::vector<score_type> matrix_storage{create_matrix_storage<score_type, matrix_order>()};

    // Note: We construct the internal data representation of the matrix depending on the matrix_major_order.
    // This will ensure that all test matrices are independent of the matrix_major_order when accessed via the same
    // matrix coordinate.
    matrix_type matrix{seqan3::detail::number_rows{3}, seqan3::detail::number_cols{4}, matrix_storage};

    template <typename value1_t, typename value2_t>
        requires !(seqan3::simd::simd_concept<std::decay_t<value1_t>> &&
                   seqan3::simd::simd_concept<std::decay_t<value2_t>>)
    static void expect_eq(value1_t v1, value2_t v2)
    {
        EXPECT_EQ(v1, v2);
    }

    template <typename value1_t, typename value2_t>
        requires (seqan3::simd::simd_concept<std::decay_t<value1_t>> &&
                  seqan3::simd::simd_concept<std::decay_t<value2_t>>)
    static void expect_eq(value1_t && v1, value2_t && v2)
    {
        SIMD_EQ(v1, v2);
    }
};

using testing_types = ::testing::Types<test_matrix_t<int, seqan3::detail::matrix_major_order::row>,
                                       test_matrix_t<int, seqan3::detail::matrix_major_order::column>,
                                       test_matrix_t<seqan3::simd::simd_type_t<int>,
                                                     seqan3::detail::matrix_major_order::row>,
                                       test_matrix_t<seqan3::simd::simd_type_t<int>,
                                                     seqan3::detail::matrix_major_order::column>>;
TYPED_TEST_SUITE(two_dimensional_matrix_test, testing_types, );

TEST(two_dimensional_matrix_test, initializer_list)
{
    [[maybe_unused]] seqan3::detail::two_dimensional_matrix<int> matrix1{seqan3::detail::number_rows{0},
                                                                         seqan3::detail::number_cols{0}, {}};
    [[maybe_unused]] seqan3::detail::two_dimensional_matrix<int> matrix2{seqan3::detail::number_rows{1},
                                                                         seqan3::detail::number_cols{1}, {0}};
}

TYPED_TEST(two_dimensional_matrix_test, concepts)
{
    using matrix_type = typename TestFixture::matrix_type;

    EXPECT_TRUE(seqan3::detail::matrix<matrix_type>);
    EXPECT_TRUE(std::ranges::input_range<matrix_type>);
    EXPECT_TRUE(std::ranges::forward_range<matrix_type>);
    EXPECT_TRUE(std::ranges::bidirectional_range<matrix_type>);
    EXPECT_TRUE(std::ranges::random_access_range<matrix_type>);
}

TYPED_TEST(two_dimensional_matrix_test, construction)
{
    using matrix_type = typename TestFixture::matrix_type;
    using score_type = typename TestFixture::score_type;

    EXPECT_TRUE(std::is_default_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_copy_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_move_constructible_v<matrix_type>);
    EXPECT_TRUE(std::is_copy_assignable_v<matrix_type>);
    EXPECT_TRUE(std::is_move_assignable_v<matrix_type>);
    EXPECT_TRUE(std::is_destructible_v<matrix_type>);

    EXPECT_TRUE((std::is_constructible_v<matrix_type,
                                         seqan3::detail::number_rows,
                                         seqan3::detail::number_cols>));
    EXPECT_TRUE((std::is_constructible_v<matrix_type,
                                         seqan3::detail::number_rows,
                                         seqan3::detail::number_cols,
                                         std::vector<score_type>>));
    EXPECT_TRUE((std::is_constructible_v<matrix_type,
                                         seqan3::detail::number_rows,
                                         seqan3::detail::number_cols,
                                         matrix_type>));
}

TYPED_TEST(two_dimensional_matrix_test, cols)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{seqan3::detail::number_rows{3}, seqan3::detail::number_cols{4}};
    EXPECT_EQ(matrix.cols(), 4u);
}

TYPED_TEST(two_dimensional_matrix_test, rows)
{
    using matrix_type = typename TestFixture::matrix_type;

    matrix_type matrix{seqan3::detail::number_rows{3}, seqan3::detail::number_cols{4}};
    EXPECT_EQ(matrix.rows(), 3u);
}

TYPED_TEST(two_dimensional_matrix_test, range)
{
    // For an explanation how this works see iterator_fixture further below in this file.
    auto it = this->matrix_storage.begin();
    for (auto cell: this->matrix)
        this->expect_eq(cell, *(it++));
}

TYPED_TEST(two_dimensional_matrix_test, subscript)
{
    // Note: Even if the internal storage has a different data layout, accessing the data via a matrix coordinate
    // yields the same cell and thus the same data.
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{0u}}],
                    this->expected_matrix_content[0]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{1u}}],
                    this->expected_matrix_content[1]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{2u}}],
                    this->expected_matrix_content[2]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{3u}}],
                    this->expected_matrix_content[3]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{0u}}],
                    this->expected_matrix_content[4]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{1u}}],
                    this->expected_matrix_content[5]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{2u}}],
                    this->expected_matrix_content[6]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{3u}}],
                    this->expected_matrix_content[7]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{2u},  seqan3::detail::column_index_type{0u}}],
                    this->expected_matrix_content[8]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{1u}}],
                    this->expected_matrix_content[9]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{2u}}],
                    this->expected_matrix_content[10]);
    this->expect_eq(this->matrix[{seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{3u}}],
                    this->expected_matrix_content[11]);
}

TYPED_TEST(two_dimensional_matrix_test, at)
{
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{0u}}),
                    this->expected_matrix_content[0]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{1u}}),
                    this->expected_matrix_content[1]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{2u}}),
                    this->expected_matrix_content[2]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{0u}, seqan3::detail::column_index_type{3u}}),
                    this->expected_matrix_content[3]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{0u}}),
                    this->expected_matrix_content[4]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{1u}}),
                    this->expected_matrix_content[5]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{2u}}),
                    this->expected_matrix_content[6]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{1u}, seqan3::detail::column_index_type{3u}}),
                    this->expected_matrix_content[7]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{0u}}),
                    this->expected_matrix_content[8]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{1u}}),
                    this->expected_matrix_content[9]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{2u}}),
                    this->expected_matrix_content[10]);
    this->expect_eq(this->matrix.at({seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{3u}}),
                    this->expected_matrix_content[11]);

    EXPECT_THROW((this->matrix.at({seqan3::detail::row_index_type{3u}, seqan3::detail::column_index_type{3u}})),
                 std::invalid_argument);
    EXPECT_THROW((this->matrix.at({seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{4u}})),
                 std::invalid_argument);
}

TYPED_TEST(two_dimensional_matrix_test, construction_other_order)
{
    // Test that changing the matrix layout works.
    static constexpr seqan3::detail::matrix_major_order other_matrix_order = TestFixture::matrix_order
                                                                           == seqan3::detail::matrix_major_order::row
                                                                           ? seqan3::detail::matrix_major_order::column
                                                                           : seqan3::detail::matrix_major_order::row;

    // Test that the implicit conversion of the underlying value_type works.
    // This does not work for every SIMD backend, in such a case we use the original score type.
    using score_type = typename TestFixture::score_type;
    using unsigned_score_type = typename make_unsigned_score_type<score_type>::type;
    using new_score_type = std::conditional_t<std::assignable_from<score_type &, unsigned_score_type &>,
                                              unsigned_score_type,
                                              score_type>;

    // The new matrix type that has a new matrix major order and might have an implicit value conversion.
    using converted_matrix_t = test_matrix_t<new_score_type, other_matrix_order>;
    converted_matrix_t converted_matrix{this->matrix};

    EXPECT_EQ(converted_matrix.rows(), this->matrix.rows());
    EXPECT_EQ(converted_matrix.cols(), this->matrix.cols());

    // Note: We changed the internal data layout, but accessing the converted matrix by a coordinate yields the same
    // cell content as accessing the original matrix by the same coordinate.
    for (unsigned row = 0; row < this->matrix.rows(); ++row)
    {
        for (unsigned col = 0; col < this->matrix.cols(); ++col)
        {
            seqan3::detail::matrix_coordinate const idx{seqan3::detail::row_index_type{row},
                                                        seqan3::detail::column_index_type{col}};
            new_score_type actual = converted_matrix[idx];
            new_score_type expected = static_cast<new_score_type>(this->matrix[idx]);
            this->expect_eq(actual, expected);
        }
    }
}

//-----------------------------------------------------------------------------
// Iterator tests
//-----------------------------------------------------------------------------

template <typename score_t, seqan3::detail::matrix_major_order order>
struct iterator_fixture<test_matrix_t<score_t, order>> : two_dimensional_matrix_test<test_matrix_t<score_t, order>>
{
    using matrix_type = test_matrix_t<score_t, order>;
    using base_t = two_dimensional_matrix_test<matrix_type>;

    // test_range is a random access range.
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    // This one is a bit tricky. The one dimensional overloads of the iterator advancing operations, like it+=5, have
    // the property that they advance the iterator in major matrix order, e.g. it+=5 in a row-wise matrix will be the
    // same as it+=matrix_offset{row{0}, col{5}} and it+=5 in a column-wise matrix will be the same as
    // it+=matrix_offset{row{5}, col{0}}. This effectively means that the iterator behaves exactly like the iterator of
    // the internal data storage, i.e. the iterator's behaviour does not depend on the major matrix order. Thus, passing
    // expected_range as data storage into test_range, will result in an equal range as test_range.
    std::vector<score_t> expected_range{this->expected_matrix_content};
    matrix_type test_range{seqan3::detail::number_rows{3}, seqan3::detail::number_cols{4}, expected_range};
};

INSTANTIATE_TYPED_TEST_SUITE_P(two_dimensional_iterator, iterator_fixture, testing_types, );

template <typename matrix_t>
struct two_dimensional_matrix_iterator_test : two_dimensional_matrix_test<matrix_t>
{
    using matrix_type = matrix_t;
    using base_t = two_dimensional_matrix_test<matrix_type>;
    using iterator_type = typename matrix_type::iterator;

    matrix_type test_range{seqan3::detail::number_rows{3},
                           seqan3::detail::number_cols{4},
                           base_t::expected_matrix_content};
};

TYPED_TEST_SUITE(two_dimensional_matrix_iterator_test, testing_types, );

TYPED_TEST(two_dimensional_matrix_iterator_test, two_dimensional_concept)
{
    EXPECT_TRUE(seqan3::detail::two_dimensional_matrix_iterator<typename TestFixture::iterator_type>);
}

TYPED_TEST(two_dimensional_matrix_iterator_test, update_by_matrix_offset_add)
{
    auto it = this->matrix.begin();
    auto it_advanced = it += seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                           seqan3::detail::column_index_type{2}};

    this->expect_eq(*it, this->expected_matrix_content[6]);
    this->expect_eq(*it_advanced, this->expected_matrix_content[6]);
}

TYPED_TEST(two_dimensional_matrix_iterator_test, advance_by_matrix_offset_add)
{
    auto it = this->matrix.begin();
    auto it_advanced = it + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                          seqan3::detail::column_index_type{2}};

    this->expect_eq(*it, this->expected_matrix_content[0]);
    this->expect_eq(*it_advanced, this->expected_matrix_content[6]);
}

TYPED_TEST(two_dimensional_matrix_iterator_test, advance_by_matrix_offset_add_friend)
{
    auto it = this->matrix.begin();
    auto it_advanced = seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                     seqan3::detail::column_index_type{2}} + it;

    this->expect_eq(*it, this->expected_matrix_content[0]);
    this->expect_eq(*it_advanced, this->expected_matrix_content[6]);
}

TYPED_TEST(two_dimensional_matrix_iterator_test, update_by_matrix_offset_subtract)
{
    auto it = this->matrix.begin() + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                   seqan3::detail::column_index_type{3}};
    auto it_advanced = it -= seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                           seqan3::detail::column_index_type{2}};

    this->expect_eq(*it, this->expected_matrix_content[5]);
    this->expect_eq(*it_advanced, this->expected_matrix_content[5]);
}

TYPED_TEST(two_dimensional_matrix_iterator_test, advance_by_matrix_offset_subtract)
{
    auto it = this->matrix.begin() + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                   seqan3::detail::column_index_type{3}};
    auto it_advanced = it - seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                          seqan3::detail::column_index_type{2}};

    this->expect_eq(*it, this->expected_matrix_content[11]);
    this->expect_eq(*it_advanced, this->expected_matrix_content[5]);
}

TYPED_TEST(two_dimensional_matrix_iterator_test, coordinate)
{
    auto it = this->matrix.begin();
    seqan3::detail::matrix_offset col_inc{seqan3::detail::row_index_type{0}, seqan3::detail::column_index_type{1}};

    // iterate row wise in the matrix
    for (std::ptrdiff_t pos = 0; it != this->matrix.end(); it += col_inc, ++pos)
    {
        auto expected_col = pos % this->matrix.cols();
        auto expected_row = pos / this->matrix.cols();
        EXPECT_EQ(it.coordinate().col, expected_col);
        EXPECT_EQ(it.coordinate().row, expected_row);
    }
}
