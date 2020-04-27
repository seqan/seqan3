// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/iterator>
#include <type_traits>
#include <vector>

#include <seqan3/alignment/matrix/detail/trace_iterator.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/range/views/to.hpp>

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

    seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> matrix{seqan3::detail::number_rows{3},
                                                                                    seqan3::detail::number_cols{4},
                                                                                    std::vector
    {
        N,           LO, L,          L,
        UO, D | LO | UO, L, D | L | UO,
        U,       LO | U, D,          L
    }};

    using trace_iterator_type = decltype(seqan3::detail::trace_iterator{matrix.begin()});
    using path_type = std::ranges::subrange<trace_iterator_type, std::ranges::default_sentinel_t>;

    path_type path(seqan3::detail::matrix_offset const & offset)
    {
        return path_type{trace_iterator_type{matrix.begin() + offset}, std::ranges::default_sentinel};
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
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                         seqan3::detail::column_index_type{3}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 5u);
    EXPECT_EQ(vec, (std::vector{L, L, L, U, U}));
}

TEST_F(trace_iterator_fixture, trace_path_2_2)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                         seqan3::detail::column_index_type{2}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 2u);
    EXPECT_EQ(vec, (std::vector{D, D}));
}

TEST_F(trace_iterator_fixture, trace_path_2_1)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                         seqan3::detail::column_index_type{1}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 3u);
    EXPECT_EQ(vec, (std::vector{U, U, L}));
}

TEST_F(trace_iterator_fixture, trace_path_2_0)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                         seqan3::detail::column_index_type{0}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 2u);
    EXPECT_EQ(vec, (std::vector{U, U}));
}

TEST_F(trace_iterator_fixture, trace_path_1_3)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                         seqan3::detail::column_index_type{3}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 3u);
    EXPECT_EQ(vec, (std::vector{D, L, L}));
}

TEST_F(trace_iterator_fixture, trace_path_1_2)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                         seqan3::detail::column_index_type{2}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 3u);
    EXPECT_EQ(vec, (std::vector{L, L, U}));
}

TEST_F(trace_iterator_fixture, trace_path_1_1)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                         seqan3::detail::column_index_type{1}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 1u);
    EXPECT_EQ(vec, (std::vector{D}));
}

TEST_F(trace_iterator_fixture, trace_path_1_0)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{1},
                                                         seqan3::detail::column_index_type{0}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 1u);
    EXPECT_EQ(vec, (std::vector{U}));
}

TEST_F(trace_iterator_fixture, trace_path_0_3)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                         seqan3::detail::column_index_type{3}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 3u);
    EXPECT_EQ(vec, (std::vector{L, L, L}));
}

TEST_F(trace_iterator_fixture, trace_path_0_2)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                         seqan3::detail::column_index_type{2}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 2u);
    EXPECT_EQ(vec, (std::vector{L, L}));
}

TEST_F(trace_iterator_fixture, trace_path_0_1)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                         seqan3::detail::column_index_type{1}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 1u);
    EXPECT_EQ(vec, (std::vector{L}));
}

TEST_F(trace_iterator_fixture, trace_path_0_0)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{0},
                                                         seqan3::detail::column_index_type{0}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 0u);
    EXPECT_EQ(vec, (std::vector<seqan3::detail::trace_directions>{}));
}

TEST_F(trace_iterator_fixture, coordinate)
{
    auto p = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                seqan3::detail::column_index_type{3}});
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
    using const_iterator_type = decltype(seqan3::detail::trace_iterator{base_t::matrix.cbegin()});

    // Test forward iterator concept.
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    struct test_range_type
    {
        auto begin()
        {
            return iterator_type{matrix.begin() + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                                seqan3::detail::column_index_type{3}}};
        }

        auto begin() const
        {
            return const_iterator_type{matrix.cbegin() +
                                       seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                     seqan3::detail::column_index_type{3}}};
        }

        auto end() { return std::ranges::default_sentinel; }
        auto end() const { return std::ranges::default_sentinel; }

        decltype(base_t::matrix) matrix;
    };

    test_range_type test_range{base_t::matrix};
    std::vector<seqan3::detail::trace_directions> expected_range{L, L, L, U, U};
};

INSTANTIATE_TYPED_TEST_SUITE_P(trace_iterator, iterator_fixture, trace_iterator_fixture, );
