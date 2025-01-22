// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iterator>
#include <type_traits>
#include <vector>

#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../../../range/iterator_test_template.hpp"

using seqan3::operator|;

struct trace_iterator_fixture : public ::testing::Test
{
    static constexpr seqan3::detail::trace_directions N = seqan3::detail::trace_directions::none;
    static constexpr seqan3::detail::trace_directions D = seqan3::detail::trace_directions::diagonal;
    static constexpr seqan3::detail::trace_directions U = seqan3::detail::trace_directions::up;
    static constexpr seqan3::detail::trace_directions UO = seqan3::detail::trace_directions::up_open;
    static constexpr seqan3::detail::trace_directions L = seqan3::detail::trace_directions::left;
    static constexpr seqan3::detail::trace_directions LO = seqan3::detail::trace_directions::left_open;

    seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> matrix{
        seqan3::detail::number_rows{3},
        seqan3::detail::number_cols{4},
        std::vector{N, LO, L, L, UO, D | LO | UO, L, D | L | UO, U, LO | U, D, L}};

    using trace_iterator_type = seqan3::detail::trace_iterator<decltype(matrix.begin())>;
    using path_type = std::ranges::subrange<trace_iterator_type, std::default_sentinel_t>;

    path_type path(seqan3::detail::matrix_offset const & offset)
    {
        return path_type{trace_iterator_type{matrix.begin() + offset}, std::default_sentinel};
    }
};

TEST_F(trace_iterator_fixture, concepts)
{
    EXPECT_TRUE(std::ranges::view<path_type>);
    EXPECT_TRUE(std::ranges::input_range<path_type>);
    EXPECT_TRUE(std::ranges::forward_range<path_type>);
    EXPECT_FALSE(std::ranges::bidirectional_range<path_type>);
}

TEST_F(trace_iterator_fixture, type_deduction)
{
    seqan3::detail::trace_iterator it{matrix.begin()};
    EXPECT_TRUE((std::is_same_v<decltype(it), seqan3::detail::trace_iterator<decltype(matrix.begin())>>));
}

TEST_F(trace_iterator_fixture, trace_path_2_3)
{
    EXPECT_RANGE_EQ(
        (std::vector{L, L, L, U, U}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2}, seqan3::detail::column_index_type{3}}));
}

TEST_F(trace_iterator_fixture, trace_path_2_2)
{
    EXPECT_RANGE_EQ(
        (std::vector{D, D}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2}, seqan3::detail::column_index_type{2}}));
}

TEST_F(trace_iterator_fixture, trace_path_2_1)
{
    EXPECT_RANGE_EQ(
        (std::vector{U, U, L}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2}, seqan3::detail::column_index_type{1}}));
}

TEST_F(trace_iterator_fixture, trace_path_2_0)
{
    EXPECT_RANGE_EQ(
        (std::vector{U, U}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2}, seqan3::detail::column_index_type{0}}));
}

TEST_F(trace_iterator_fixture, trace_path_1_3)
{
    EXPECT_RANGE_EQ(
        (std::vector{D, L, L}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1}, seqan3::detail::column_index_type{3}}));
}

TEST_F(trace_iterator_fixture, trace_path_1_2)
{
    EXPECT_RANGE_EQ(
        (std::vector{L, L, U}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1}, seqan3::detail::column_index_type{2}}));
}

TEST_F(trace_iterator_fixture, trace_path_1_1)
{
    EXPECT_RANGE_EQ(
        (std::vector{D}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1}, seqan3::detail::column_index_type{1}}));
}

TEST_F(trace_iterator_fixture, trace_path_1_0)
{
    EXPECT_RANGE_EQ(
        (std::vector{U}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1}, seqan3::detail::column_index_type{0}}));
}

TEST_F(trace_iterator_fixture, trace_path_0_3)
{
    EXPECT_RANGE_EQ(
        (std::vector{L, L, L}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0}, seqan3::detail::column_index_type{3}}));
}

TEST_F(trace_iterator_fixture, trace_path_0_2)
{
    EXPECT_RANGE_EQ(
        (std::vector{L, L}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0}, seqan3::detail::column_index_type{2}}));
}

TEST_F(trace_iterator_fixture, trace_path_0_1)
{
    EXPECT_RANGE_EQ(
        (std::vector{L}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0}, seqan3::detail::column_index_type{1}}));
}

TEST_F(trace_iterator_fixture, trace_path_0_0)
{
    EXPECT_RANGE_EQ(
        (std::vector<seqan3::detail::trace_directions>{}),
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0}, seqan3::detail::column_index_type{0}}));
}

TEST_F(trace_iterator_fixture, coordinate)
{
    auto p =
        path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2}, seqan3::detail::column_index_type{3}});
    auto it = p.begin();

    EXPECT_EQ(it.coordinate().row, 2u);
    EXPECT_EQ(it.coordinate().col, 3u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 2u);
    EXPECT_EQ(it.coordinate().col, 2u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 2u);
    EXPECT_EQ(it.coordinate().col, 1u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 2u);
    EXPECT_EQ(it.coordinate().col, 0u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 1u);
    EXPECT_EQ(it.coordinate().col, 0u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 0u);
    EXPECT_EQ(it.coordinate().col, 0u);
}

//-----------------------------------------------------------------------------
// Iterator tests
//-----------------------------------------------------------------------------

template <>
struct iterator_fixture<trace_iterator_fixture> : public trace_iterator_fixture
{
    using base_t = trace_iterator_fixture;

    using iterator_type = typename base_t::trace_iterator_type;
    using const_iterator_type = seqan3::detail::trace_iterator<decltype(base_t::matrix.cbegin())>;

    // Test forward iterator concept.
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    struct test_range_type
    {
        auto begin()
        {
            return iterator_type{matrix.begin()
                                 + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                 seqan3::detail::column_index_type{3}}};
        }

        auto begin() const
        {
            return const_iterator_type{matrix.cbegin()
                                       + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                       seqan3::detail::column_index_type{3}}};
        }

        auto end()
        {
            return std::default_sentinel;
        }
        auto end() const
        {
            return std::default_sentinel;
        }

        decltype(base_t::matrix) matrix;
    };

    test_range_type test_range{base_t::matrix};
    std::vector<seqan3::detail::trace_directions> expected_range{L, L, L, U, U};
};

INSTANTIATE_TYPED_TEST_SUITE_P(trace_iterator, iterator_fixture, trace_iterator_fixture, );
