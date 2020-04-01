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

#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/search/search_result_range.hpp>

#include "../range/iterator_test_template.hpp"

// ----------------------------------------------------------------------------
// Simple executor used as mock for the test.
// ----------------------------------------------------------------------------

struct dummy_search_algorithm
{
    using value_type      = size_t;
    using reference       = value_type;
    using difference_type = std::ptrdiff_t;

    template <typename query_t>
    auto operator()(query_t && /*query*/)
    {
        return std::vector<size_t>{0, 1, 2, 3, 4};
    }
};

// ----------------------------------------------------------------------------
// Testing iterator.
// ----------------------------------------------------------------------------

using search_result_range_t = seqan3::search_result_range<dummy_search_algorithm,
                                                          seqan3::type_reduce_view<std::vector<std::string>& >>;
using search_range_iterator = std::ranges::iterator_t<search_result_range_t>;

struct search_result_range_test : ::testing::Test
{
    std::vector<std::string> query_range{std::string{"query1"}, std::string{"query1"}, std::string{"query3"}};
};

template <>
struct iterator_fixture<search_range_iterator> : search_result_range_test
{
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    search_result_range_t test_range{dummy_search_algorithm{}, this->query_range | seqan3::views::type_reduce};
    std::vector<std::pair<size_t, size_t>> expected_range{{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
                                                          {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
                                                          {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}};
};

INSTANTIATE_TYPED_TEST_SUITE_P(search_range_iterator, iterator_fixture, search_range_iterator, );

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

    EXPECT_TRUE((std::is_constructible_v<search_result_range_t,
                                         dummy_search_algorithm,
                                         seqan3::type_reduce_view<std::vector<std::string>& >>));
}

TEST_F(search_result_range_test, type_deduction)
{
    seqan3::search_result_range rng{dummy_search_algorithm{}, query_range | seqan3::views::type_reduce};
    EXPECT_TRUE((std::is_same_v<decltype(rng), search_result_range_t>));
}

TEST_F(search_result_range_test, empty_query_range)
{
    query_range.clear();
    seqan3::search_result_range rng{dummy_search_algorithm{}, query_range | seqan3::views::type_reduce};

    EXPECT_TRUE(rng.begin() == rng.end());
}
