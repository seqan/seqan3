// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/matrix/detail/alignment_matrix_column_major_range_base.hpp>

class test_matrix : public seqan3::detail::alignment_matrix_column_major_range_base<test_matrix>
{
public:

    using base_t = seqan3::detail::alignment_matrix_column_major_range_base<test_matrix>;

    using element_type = int;
    using alignment_column_type = typename base_t::alignment_column_type;
    using column_data_view_type = std::span<element_type>;

    struct proxy_type
    {
        element_type & value;
    };

    using value_type = proxy_type;
    using reference = proxy_type;

    test_matrix() = default;
    test_matrix(test_matrix const &) = default;
    test_matrix(test_matrix &&) = default;
    test_matrix & operator=(test_matrix const &) = default;
    test_matrix & operator=(test_matrix &&) = default;
    ~test_matrix() = default;

    std::vector<element_type> data{};
    size_t num_cols{};
    size_t num_rows{};
    size_t num_create{};
    size_t num_pre{};
    size_t num_post{};

    alignment_column_type initialise_column(size_t const column_index) noexcept
    {
        return alignment_column_type{*this,
                                     column_data_view_type{std::addressof(data[num_rows * column_index]), num_rows}};
    }

    template <std::random_access_iterator iter_t>
    constexpr reference make_proxy(iter_t iter) noexcept
    {
        return {*iter};
    }

    template <std::random_access_iterator iter_t>
    constexpr void on_column_iterator_creation(iter_t) noexcept
    {
        ++num_create;
    }

    template <std::random_access_iterator iter_t>
    constexpr void before_column_iterator_increment(iter_t ) noexcept
    {
        ++num_pre;
    }

    template <std::random_access_iterator iter_t>
    constexpr void after_column_iterator_increment(iter_t ) noexcept
    {
        ++num_post;
    }
};

class alignment_matrix_column_major_range_base_test : public ::testing::Test
{
protected:

    alignment_matrix_column_major_range_base_test()
    {
        matrix.num_cols = 4;
        matrix.num_rows = 5;
        matrix.data.resize(matrix.num_cols * matrix.num_rows);

        int k = 0;
        for (unsigned i = 0; i < matrix.num_cols; ++i)
            for (unsigned j = 0; j < matrix.num_rows; ++j, ++k)
                matrix.data[i * matrix.num_rows + j] = k;
    }

    test_matrix matrix{};
};

TEST(alignment_matrix_column_major_range_base, concepts)
{
    using outer_it = std::ranges::iterator_t<test_matrix>;
    using column_t = std::iter_value_t<outer_it>;
    using inner_it = std::ranges::iterator_t<column_t>;

    EXPECT_TRUE(std::ranges::input_range<test_matrix>);
    EXPECT_TRUE(std::input_iterator<outer_it>);
    EXPECT_TRUE(std::input_iterator<inner_it>);
    EXPECT_TRUE(std::ranges::input_range<column_t>);
    EXPECT_TRUE(std::ranges::view<column_t>);
}

TEST_F(alignment_matrix_column_major_range_base_test, begin_end)
{
    auto it = matrix.begin();
    auto col = *it;
    auto col_it = col.begin();
    EXPECT_EQ((*col_it).value, 0);
    EXPECT_NE(col_it, col.end());
    EXPECT_NE(it, matrix.end());
}

TEST_F(alignment_matrix_column_major_range_base_test, iterate_columns)
{
    size_t count = 0;
    for (auto it = matrix.begin(); it != matrix.end(); ++it, ++count);
    EXPECT_EQ(count, 4u);
}

TEST_F(alignment_matrix_column_major_range_base_test, iterate_num_rows)
{
    auto it = matrix.begin();
    auto col = *it;
    size_t count = 0;
    for (auto it = col.begin(); it != col.end(); ++it, ++count);
    EXPECT_EQ(count, 5u);
}

TEST_F(alignment_matrix_column_major_range_base_test, iterate_matrix)
{
    auto mat_it = matrix.begin();
    auto col = *mat_it;
    int cmp = 0;
    for (auto it = col.begin(); it != col.end(); ++it, ++cmp)
        EXPECT_EQ((*it).value, cmp);

    ++mat_it;
    col = *mat_it;
    for (auto it = col.begin(); it != col.end(); ++it, ++cmp)
        EXPECT_EQ((*it).value, cmp);

    ++mat_it;
    col = *mat_it;
    for (auto it = col.begin(); it != col.end(); ++it, ++cmp)
        EXPECT_EQ((*it).value, cmp);

    ++mat_it;
    col = *mat_it;
    for (auto it = col.begin(); it != col.end(); ++it, ++cmp)
        EXPECT_EQ((*it).value, cmp);

    EXPECT_EQ(matrix.num_create, 4u);
    EXPECT_EQ(matrix.num_pre, 20u);
    EXPECT_EQ(matrix.num_post, 20u);
}
