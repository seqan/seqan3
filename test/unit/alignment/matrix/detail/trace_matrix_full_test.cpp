// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/alignment/matrix/detail/trace_matrix_full.hpp>

#include "../../../range/iterator_test_template.hpp"

using trace_t = seqan3::detail::trace_directions;
using matrix_t = seqan3::detail::trace_matrix_full<trace_t>;
using matrix_iterator_t = std::ranges::iterator_t<matrix_t>;

template <>
struct iterator_fixture<matrix_iterator_t> : public ::testing::Test
{
    using materialised_column_t = std::vector<std::tuple<trace_t, trace_t, trace_t>>;

    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    static constexpr seqan3::detail::trace_directions none = seqan3::detail::trace_directions::none;

    // Single column with 5 entries as the seq2 has size 4 (need one more for the initialisation row).
    materialised_column_t column = materialised_column_t{5, std::tuple{none, none, none}};
    std::vector<materialised_column_t> expected_range{column, column, column, column};
    matrix_t test_range;

    void SetUp()
    {
        std::string seq1 = "abc";
        std::string seq2 = "abcd";

        test_range.resize(seqan3::detail::column_index_type<size_t>{4}, seqan3::detail::row_index_type<size_t>{5});
    }

    template <typename actual_column_t, typename expected_column_t>
    static void expect_eq(actual_column_t && actual_column, expected_column_t && expected_column)
    {
        auto actual_it = actual_column.begin();
        auto expected_it = expected_column.begin();
        for (; actual_it != actual_column.end(); ++actual_it, ++expected_it)
        {
            using std::get;
            auto actual_cell = *actual_it;
            auto expected_cell = *expected_it;

            EXPECT_EQ(get<0>(actual_cell), get<0>(expected_cell));
            EXPECT_EQ(get<1>(actual_cell), get<1>(expected_cell));
            EXPECT_EQ(get<2>(actual_cell), get<2>(expected_cell));
        }
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(trace_matrix_full_test, iterator_fixture, matrix_iterator_t, );

TEST(trace_matrix_full_test, viewable_range_proxy)
{
    EXPECT_TRUE(std::ranges::view<std::iter_reference_t<matrix_iterator_t>>);
}

TEST(trace_matrix_full_test, trace_path)
{
    matrix_t matrix{};
    matrix.resize(seqan3::detail::column_index_type<size_t>{4}, seqan3::detail::row_index_type<size_t>{3});
    auto trace_column_it = matrix.begin();
    auto trace_column = *trace_column_it;

    // Initialise column 0
    auto trace_cell_it = trace_column.begin();
    *trace_cell_it = std::tuple{trace_t::none, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::up_open, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::up, trace_t::none, trace_t::none};

    // Initialise column 1
    trace_column = *++trace_column_it;
    trace_cell_it = trace_column.begin();
    *trace_cell_it = std::tuple{trace_t::left_open, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::diagonal, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::up_open, trace_t::none, trace_t::none};

    // Initialise column 2
    trace_column = *++trace_column_it;
    trace_cell_it = trace_column.begin();
    *trace_cell_it = std::tuple{trace_t::left, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::diagonal, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::left_open, trace_t::none, trace_t::none};

    // Initialise column 3
    trace_column = *++trace_column_it;
    trace_cell_it = trace_column.begin();
    *trace_cell_it = std::tuple{trace_t::left, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::up_open, trace_t::none, trace_t::none};
    *++trace_cell_it = std::tuple{trace_t::left, trace_t::none, trace_t::none};

    EXPECT_TRUE(++trace_cell_it == trace_column.end());
    EXPECT_TRUE(++trace_column_it == matrix.end());

    auto trace_path = matrix.trace_path(
        seqan3::detail::matrix_coordinate{seqan3::detail::row_index_type{2u}, seqan3::detail::column_index_type{3u}});

    auto trace_path_it = trace_path.begin();
    EXPECT_EQ(*trace_path_it, trace_t::left);
    EXPECT_EQ(*++trace_path_it, trace_t::left);
    EXPECT_EQ(*++trace_path_it, trace_t::up);
    EXPECT_EQ(*++trace_path_it, trace_t::diagonal);
    EXPECT_EQ(*++trace_path_it, trace_t::none);
    EXPECT_TRUE(trace_path_it == trace_path.end());
}

TEST(trace_matrix_full_test, invalid_trace_path_coordinate)
{
    matrix_t matrix{};
    matrix.resize(seqan3::detail::column_index_type<size_t>{4}, seqan3::detail::row_index_type<size_t>{3});

    EXPECT_THROW((matrix.trace_path(seqan3::detail::matrix_coordinate{seqan3::detail::row_index_type{3u},
                                                                      seqan3::detail::column_index_type{3u}})),
                 std::invalid_argument);
    EXPECT_THROW((matrix.trace_path(seqan3::detail::matrix_coordinate{seqan3::detail::row_index_type{2u},
                                                                      seqan3::detail::column_index_type{4u}})),
                 std::invalid_argument);
}
