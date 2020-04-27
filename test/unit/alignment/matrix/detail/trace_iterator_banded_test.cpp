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

#include <seqan3/alignment/matrix/detail/trace_iterator_banded.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/range/views/to.hpp>

#include "../../../range/iterator_test_template.hpp"

struct trace_iterator_banded_test : public ::testing::Test
{
    static constexpr seqan3::detail::trace_directions N = seqan3::detail::trace_directions::none;
    static constexpr seqan3::detail::trace_directions D = seqan3::detail::trace_directions::diagonal;
    static constexpr seqan3::detail::trace_directions U = seqan3::detail::trace_directions::up;
    static constexpr seqan3::detail::trace_directions UO = seqan3::detail::trace_directions::up_open;
    static constexpr seqan3::detail::trace_directions L = seqan3::detail::trace_directions::left;
    static constexpr seqan3::detail::trace_directions LO = seqan3::detail::trace_directions::left_open;

    seqan3::detail::two_dimensional_matrix<seqan3::detail::trace_directions> matrix{seqan3::detail::number_rows{5},
                                                                                    seqan3::detail::number_cols{6},
                                                                                    std::vector
    {  // emulating band {row:2, col:2}
         N,  N,  L,  D,  D,  L,
         N, LO,  D,  D, UO, UO,
         N,  D,  D, LO,  D,  U,
        UO,  D,  D,  D,  D,  N,
         U,  D,  D,  D,  N,  N
    }};
            // Band view
        //   0  1  2  3  4  5
        //0  N LO  L
        //1 UO  D  D  D
        //2  U  D  D  D  D
        //3     D  D LO UO  L
        //4        D  D  D UO
        //5           D  D  U

    using trace_iterator_type = decltype(seqan3::detail::trace_iterator_banded{matrix.begin(),
                                                                               seqan3::detail::column_index_type{0}});
    using path_type = std::ranges::subrange<trace_iterator_type, std::ranges::default_sentinel_t>;

    path_type path(seqan3::detail::matrix_offset const & offset)
    {
        return path_type{trace_iterator_type{matrix.begin() + offset, seqan3::detail::column_index_type{2}},
                         std::ranges::default_sentinel};
    }
};

TEST_F(trace_iterator_banded_test, concepts)
{
    EXPECT_TRUE(std::ranges::view<path_type>);
    EXPECT_TRUE(std::ranges::input_range<path_type>);
    EXPECT_TRUE(std::ranges::forward_range<path_type>);
    EXPECT_FALSE(std::ranges::bidirectional_range<path_type>);
}

TEST_F(trace_iterator_banded_test, type_deduction)
{
    seqan3::detail::trace_iterator_banded it{matrix.begin(), seqan3::detail::column_index_type{0u}};
    EXPECT_TRUE((std::is_same_v<decltype(it), seqan3::detail::trace_iterator_banded<decltype(matrix.begin())>>));
}

TEST_F(trace_iterator_banded_test, trace_path_2_5)
{
    std::vector vec = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                         seqan3::detail::column_index_type{5}})
                    | seqan3::views::to<std::vector>;

    EXPECT_EQ(vec.size(), 8u);
    EXPECT_EQ(vec, (std::vector{U, U, L, L, L, D, D, U}));
}

TEST_F(trace_iterator_banded_test, coordinate)
{
    auto p = path(seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                seqan3::detail::column_index_type{5}});
    auto it = p.begin();

    EXPECT_EQ(it.coordinate().row, 5u);
    EXPECT_EQ(it.coordinate().col, 5u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 4u);
    EXPECT_EQ(it.coordinate().col, 5u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 3u);
    EXPECT_EQ(it.coordinate().col, 5u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 3u);
    EXPECT_EQ(it.coordinate().col, 4u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 3u);
    EXPECT_EQ(it.coordinate().col, 3u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 3u);
    EXPECT_EQ(it.coordinate().col, 2u);
    ++it;
    EXPECT_EQ(it.coordinate().row, 2u);
    EXPECT_EQ(it.coordinate().col, 1u);
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
struct iterator_fixture<trace_iterator_banded_test> : public trace_iterator_banded_test
{
    using base_t = trace_iterator_banded_test;

    using iterator_type = typename base_t::trace_iterator_type;
    using const_iterator_type = decltype(seqan3::detail::trace_iterator_banded{base_t::matrix.cbegin(),
                                                                               seqan3::detail::column_index_type{0}});

    // Test forward iterator concept.
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    struct test_range_type
    {
        auto begin()
        {
            return iterator_type{matrix.begin() + seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                                seqan3::detail::column_index_type{5}},
                                 seqan3::detail::column_index_type{2}};
        }

        auto begin() const
        {
            return const_iterator_type{matrix.cbegin() +
                                       seqan3::detail::matrix_offset{seqan3::detail::row_index_type{2},
                                                                     seqan3::detail::column_index_type{5}},
                                       seqan3::detail::column_index_type{2}};
        }

        auto end() { return std::ranges::default_sentinel; }
        auto end() const { return std::ranges::default_sentinel; }

        decltype(base_t::matrix) matrix;
    };

    test_range_type test_range{base_t::matrix};
    std::vector<seqan3::detail::trace_directions> expected_range{U, U, L, L, L, D, D, U};
};

INSTANTIATE_TYPED_TEST_SUITE_P(trace_iterator_banded, iterator_fixture, trace_iterator_banded_test, );
