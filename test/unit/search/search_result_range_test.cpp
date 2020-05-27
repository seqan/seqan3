// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <optional>
#include <string>
#include <vector>

#include <seqan3/search/search_result_range.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../range/iterator_test_template.hpp"

// ----------------------------------------------------------------------------
// Simple executor used as mock for the test.
// ----------------------------------------------------------------------------

using dummy_result_type = std::pair<size_t, size_t>;

struct dummy_search_algorithm
{
    using value_type      = size_t;
    using reference       = value_type;
    using difference_type = std::ptrdiff_t;

    template <typename indexed_query_t, typename callback_t>
    auto operator()(indexed_query_t && indexed_query, callback_t && callback)
    {
        auto idx = std::get<0>(indexed_query);
        for (size_t i : std::vector<size_t>{0, 1, 2, 3, 4})
            callback(dummy_result_type{idx, i});
    }
};

// ----------------------------------------------------------------------------
// Testing iterator.
// ----------------------------------------------------------------------------

struct search_result_range_test : ::testing::Test
{
    using indexed_queries_t = std::vector<std::pair<size_t, std::string>>;
    using algorithm_result_t = std::pair<size_t, size_t>;
    using executor_t = seqan3::detail::algorithm_executor_blocking<indexed_queries_t &,
                                                                   dummy_search_algorithm,
                                                                   algorithm_result_t>;
    using search_result_range_t = seqan3::search_result_range<executor_t>;

    indexed_queries_t indexed_queries{{0, std::string{"query1"}},
                                      {1, std::string{"query1"}},
                                      {2, std::string{"query3"}}};

};

template <>
struct iterator_fixture<search_result_range_test> : search_result_range_test
{
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    search_result_range_t test_range{executor_t{this->indexed_queries, dummy_search_algorithm{}, algorithm_result_t{}}};
    std::vector<dummy_result_type> expected_range{{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
                                                  {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
                                                  {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}};
};

INSTANTIATE_TYPED_TEST_SUITE_P(search_range_iterator, iterator_fixture, search_result_range_test, );

// ----------------------------------------------------------------------------
// Testing alignment range concepts and interfaces.
// ----------------------------------------------------------------------------

TEST_F(search_result_range_test, concept_test)
{
    EXPECT_TRUE(std::ranges::input_range<search_result_range_t>);
    EXPECT_FALSE(std::ranges::forward_range<search_result_range_t>);
}

TEST_F(search_result_range_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<search_result_range_t>);
    EXPECT_FALSE(std::is_copy_constructible_v<search_result_range_t>);
    EXPECT_TRUE(std::is_move_constructible_v<search_result_range_t>);
    EXPECT_FALSE(std::is_copy_assignable_v<search_result_range_t>);
    EXPECT_TRUE(std::is_move_assignable_v<search_result_range_t>);

    EXPECT_TRUE((std::is_constructible_v<search_result_range_t, executor_t>));
}

TEST_F(search_result_range_test, type_deduction)
{
    seqan3::search_result_range rng{executor_t{indexed_queries, dummy_search_algorithm{}, algorithm_result_t{}}};
    EXPECT_TRUE((std::is_same_v<decltype(rng), search_result_range_t>));
}

TEST_F(search_result_range_test, empty_query_range)
{
    indexed_queries.clear();
    seqan3::search_result_range rng{executor_t{indexed_queries, dummy_search_algorithm{}, algorithm_result_t{}}};

    EXPECT_TRUE(rng.begin() == rng.end());
}

TEST_F(search_result_range_test, issue1799)
{
    std::vector<dummy_result_type> expected_range{{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
                                                  {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
                                                  {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}};
    { // move construction
        seqan3::search_result_range rng{executor_t{indexed_queries, dummy_search_algorithm{}, algorithm_result_t{}}};
        seqan3::search_result_range moved_range{std::move(rng)};

        EXPECT_RANGE_EQ(expected_range, moved_range);
    }

    { // move assignment
        seqan3::search_result_range rng{executor_t{indexed_queries, dummy_search_algorithm{}, algorithm_result_t{}}};
        auto && moved_range = std::move(rng);

        EXPECT_RANGE_EQ(expected_range, moved_range);
    }
}
