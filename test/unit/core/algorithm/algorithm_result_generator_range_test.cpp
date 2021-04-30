// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <optional>
#include <vector>

#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

#include "../../range/iterator_test_template.hpp"

// ----------------------------------------------------------------------------
// Simple executor used as mock for the test.
// ----------------------------------------------------------------------------

struct dummy_executor
{
    using value_type      = size_t;
    using reference       = value_type;
    using difference_type = std::ptrdiff_t;

    std::optional<size_t> next_result()
    {
        auto it = std::ranges::begin(generator);
        if (it == std::ranges::end(generator))
        {
            return {};
        }
        else
        {
            std::optional<size_t> opt{*it};
            ++it;
            return opt;
        }
    }

private:

    seqan3::detail::single_pass_input_view<decltype(std::views::iota(0u, 10u))> generator{std::views::iota(0u, 10u)};
};

// ----------------------------------------------------------------------------
// Testing iterator.
// ----------------------------------------------------------------------------

using algorithm_result_generator_range_t = seqan3::algorithm_result_generator_range<dummy_executor>;
using algorithm_result_generator_range_iterator = std::ranges::iterator_t<algorithm_result_generator_range_t>;

template <>
struct iterator_fixture<algorithm_result_generator_range_iterator> : ::testing::Test
{
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    algorithm_result_generator_range_t test_range{dummy_executor{}};
    std::vector<size_t> expected_range{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
};

INSTANTIATE_TYPED_TEST_SUITE_P(algorithm_result_generator_range_iterator, iterator_fixture, algorithm_result_generator_range_iterator, );

// ----------------------------------------------------------------------------
// Testing alignment range concepts and interfaces.
// ----------------------------------------------------------------------------

TEST(algorithm_result_generator_range, concept_test)
{
    EXPECT_TRUE(std::ranges::input_range<algorithm_result_generator_range_t>);
    EXPECT_FALSE(std::ranges::forward_range<algorithm_result_generator_range_t>);
}

TEST(algorithm_result_generator_range, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<algorithm_result_generator_range_t>);
    EXPECT_FALSE(std::is_copy_constructible_v<algorithm_result_generator_range_t>);
    EXPECT_TRUE(std::is_move_constructible_v<algorithm_result_generator_range_t>);
    EXPECT_FALSE(std::is_copy_assignable_v<algorithm_result_generator_range_t>);
    EXPECT_TRUE(std::is_move_assignable_v<algorithm_result_generator_range_t>);

    EXPECT_TRUE((std::is_constructible_v<algorithm_result_generator_range_t, dummy_executor>));
}

TEST(algorithm_result_generator_range, type_deduction)
{
    seqan3::algorithm_result_generator_range rng{dummy_executor{}};
    EXPECT_TRUE((std::is_same_v<decltype(rng), seqan3::algorithm_result_generator_range<dummy_executor>>));
}

TEST(algorithm_result_generator_range, begin)
{
    seqan3::algorithm_result_generator_range rng{dummy_executor{}};
    auto it = rng.begin();
    EXPECT_EQ(*it, 0u);
}

TEST(algorithm_result_generator_range, end)
{
    seqan3::algorithm_result_generator_range rng{dummy_executor{}};
    auto it = rng.end();
    EXPECT_FALSE(it == rng.begin());
    EXPECT_FALSE(rng.begin() == it);
}

TEST(algorithm_result_generator_range, iterable)
{
    seqan3::algorithm_result_generator_range rng{dummy_executor{}};
    size_t sum = 0;
    for (size_t res : rng)
        sum += res;

    EXPECT_EQ(sum, 45u);
}

TEST(algorithm_result_generator_range, default_construction)
{
    seqan3::algorithm_result_generator_range<dummy_executor> rng{};
    EXPECT_THROW(rng.begin(), std::runtime_error);
}
