// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>
#include <utility>

#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>

#include "../../../range/iterator_test_template.hpp"
#include "alignment_matrix_base_test_template.hpp"

using trace_matrix_t =
    std::pair<seqan3::detail::alignment_trace_matrix_full<seqan3::detail::trace_directions>, std::false_type>;
using coo_matrix_t =
    std::pair<seqan3::detail::alignment_trace_matrix_full<seqan3::detail::trace_directions, true>, std::false_type>;

using testing_types = ::testing::Types<trace_matrix_t, coo_matrix_t>;

INSTANTIATE_TYPED_TEST_SUITE_P(full_matrix, alignment_matrix_base_test, testing_types, );

//-----------------------------------------------------------------------------
// Test outer iterator
//-----------------------------------------------------------------------------
template <typename test_type>
struct outer_iterator
{};

template <typename test_type>
struct iterator_fixture<outer_iterator<test_type>> : alignment_matrix_base_test<test_type>
{
    static constexpr seqan3::detail::trace_directions N = seqan3::detail::trace_directions::none;

    using base_t = alignment_matrix_base_test<test_type>;
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    using base_t::test_range;
    std::vector<std::pair<std::pair<size_t, size_t>, seqan3::detail::trace_directions>> expected_range{{{0, 0}, N},
                                                                                                       {{0, 1}, N},
                                                                                                       {{0, 2}, N},
                                                                                                       {{0, 3}, N},
                                                                                                       {{0, 4}, N}};

    template <typename lhs_t, typename rhs_t>
    static void test(lhs_t lhs, rhs_t rhs)
    {
        EXPECT_EQ(lhs.coordinate.second, std::get<0>(std::get<0>(rhs)));
        EXPECT_EQ(lhs.coordinate.first, std::get<1>(std::get<0>(rhs)));

        if constexpr (std::is_same_v<std::tuple_element_t<0, test_type>,
                                     seqan3::detail::alignment_trace_matrix_full<seqan3::detail::trace_directions>>)
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

INSTANTIATE_TYPED_TEST_SUITE_P(trace_matrix_outer_iterator, iterator_fixture, testing_types_outer, );

//-----------------------------------------------------------------------------
// Test inner iterator
//-----------------------------------------------------------------------------
template <typename test_type>
struct inner_iterator
{};

template <typename test_type>
struct iterator_fixture<inner_iterator<test_type>> : iterator_fixture<outer_iterator<test_type>>
{
    using base_t = iterator_fixture<outer_iterator<test_type>>;
    using base_t::N;

    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    decltype(*base_t::test_range.begin()) test_range = *base_t::test_range.begin();
    std::vector<std::pair<std::pair<size_t, size_t>, seqan3::detail::trace_directions>> expected_range{{{0, 0}, N},
                                                                                                       {{1, 0}, N},
                                                                                                       {{2, 0}, N},
                                                                                                       {{3, 0}, N},
                                                                                                       {{4, 0}, N}};

    template <typename lhs_t, typename rhs_t>
    static void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        base_t::test(lhs, rhs);
    }
};

using testing_types_inner = ::testing::Types<inner_iterator<trace_matrix_t>, inner_iterator<coo_matrix_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(trace_matrix_inner_iterator, iterator_fixture, testing_types_inner, );

TEST(trace_matrix, trace_path)
{
    seqan3::detail::alignment_trace_matrix_full<seqan3::detail::trace_directions> matrix{"acgt", "acgt"};

    EXPECT_THROW((matrix.trace_path(seqan3::detail::matrix_coordinate{seqan3::detail::row_index_type{6u},
                                                                      seqan3::detail::column_index_type{4u}})),
                 std::invalid_argument);

    EXPECT_THROW((matrix.trace_path(seqan3::detail::matrix_coordinate{seqan3::detail::row_index_type{4u},
                                                                      seqan3::detail::column_index_type{6u}})),
                 std::invalid_argument);

    auto path = matrix.trace_path(
        seqan3::detail::matrix_coordinate{seqan3::detail::row_index_type{4u}, seqan3::detail::column_index_type{4u}});

    EXPECT_TRUE(path.empty());
}
