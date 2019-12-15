// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>
#include <utility>

#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full_banded.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>

#include "alignment_matrix_base_test_template.hpp"
#include "../../../range/iterator_test_template.hpp"

using namespace seqan3;

using trace_matrix_t = std::pair<detail::alignment_trace_matrix_full_banded<trace_directions>, std::true_type>;
using coo_matrix_t = std::pair<detail::alignment_trace_matrix_full_banded<trace_directions, true>, std::true_type>;

using testing_types = ::testing::Types<trace_matrix_t, coo_matrix_t>;
INSTANTIATE_TYPED_TEST_SUITE_P(full_matrix_banded,
                               alignment_matrix_base_test,
                               testing_types, );

//-----------------------------------------------------------------------------
// Test outer iterator
//-----------------------------------------------------------------------------
template <typename test_type>
struct outer_iterator
{};

template <typename test_type>
struct iterator_fixture<outer_iterator<test_type>> : alignment_matrix_base_test<test_type>
{
    static constexpr trace_directions N = trace_directions::none;

    using base_t = alignment_matrix_base_test<test_type>;
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    using base_t::test_range;
    std::vector<std::pair<std::pair<size_t, size_t>, trace_directions>> expected_range{{{2, 0}, N},
                                                                                       {{1, 1}, N},
                                                                                       {{0, 2}, N},
                                                                                       {{0, 3}, N},
                                                                                       {{0, 4}, N}};

    template <typename lhs_t, typename rhs_t>
    static void test(lhs_t lhs, rhs_t rhs)
    {
        EXPECT_EQ(lhs.coordinate.second, std::get<0>(std::get<0>(rhs)));
        EXPECT_EQ(lhs.coordinate.first, std::get<1>(std::get<0>(rhs)));

        if constexpr (std::is_same_v<std::tuple_element_t<0, test_type>,
                                     detail::alignment_trace_matrix_full_banded<trace_directions>>)
        {
            EXPECT_EQ(lhs.current, std::get<1>(rhs));
        }
    }

    template <typename lhs_t, typename rhs_t>
    static void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        test(*lhs.begin(), rhs);
    }
};

using testing_types_outer = ::testing::Types<outer_iterator<trace_matrix_t>, outer_iterator<coo_matrix_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(banded_trace_matrix_outer_iterator,
                               iterator_fixture,
                               testing_types_outer, );

//-----------------------------------------------------------------------------
// Test inner iterator
//-----------------------------------------------------------------------------
template <typename test_type>
struct inner_iterator
{};

template <typename test_type>
struct iterator_fixture<inner_iterator<test_type>> : iterator_fixture<outer_iterator<test_type>>
{
    static constexpr trace_directions N = trace_directions::none;

    using base_t = iterator_fixture<outer_iterator<test_type>>;

    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    decltype(*base_t::test_range.begin()) test_range = *base_t::test_range.begin();
    std::vector<std::pair<std::pair<size_t, size_t>, trace_directions>> expected_range{{{2, 0}, base_t::N},
                                                                                       {{3, 0}, base_t::N},
                                                                                       {{4, 0}, base_t::N}};

    template <typename lhs_t, typename rhs_t>
    static void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        base_t::test(lhs, rhs);
    }
};

using testing_types_inner = ::testing::Types<inner_iterator<trace_matrix_t>, inner_iterator<coo_matrix_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(banded_trace_matrix_inner_iterator,
                               iterator_fixture,
                               testing_types_inner, );

TEST(trace_matrix, trace_path)
{
    detail::alignment_trace_matrix_full_banded<trace_directions>
        matrix{"acgt", "acgt", static_band{lower_bound{-3}, upper_bound{3}}};

    EXPECT_THROW((matrix.trace_path(matrix_coordinate{row_index_type{7u}, column_index_type{4u}})),
                 std::invalid_argument);

    EXPECT_THROW((matrix.trace_path(matrix_coordinate{row_index_type{4u}, column_index_type{7u}})),
                 std::invalid_argument);

    auto path = matrix.trace_path(matrix_coordinate{row_index_type{4u}, column_index_type{4u}});

    EXPECT_TRUE(path.empty());
}
