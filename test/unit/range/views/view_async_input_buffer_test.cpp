// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/expect_range_eq.hpp>

#include "../iterator_test_template.hpp"

using seqan3::operator""_dna4;

using iterator_type = std::ranges::iterator_t<
    decltype(std::declval<seqan3::dna4_vector&>() | seqan3::views::async_input_buffer(3))>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::input_iterator_tag;
    static constexpr bool const_iterable = false;

    seqan3::dna4_vector vec{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};
    seqan3::dna4_vector expected_range{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};
    decltype(seqan3::views::async_input_buffer(vec, 3)) test_range = seqan3::views::async_input_buffer(vec, 3);
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

TEST(async_input_buffer, in_out)
{
    seqan3::dna4_vector vec{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};

    auto v = vec | seqan3::views::async_input_buffer(3);

    EXPECT_RANGE_EQ(vec, v);
}

TEST(async_input_buffer, in_out_empty)
{
    seqan3::dna4_vector vec{};

    auto v = vec | seqan3::views::async_input_buffer(3);

    EXPECT_TRUE(v.begin() == v.end());
}

TEST(async_input_buffer, buffer_size_zero)
{
    seqan3::dna4_vector vec{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};

    EXPECT_THROW(vec | seqan3::views::async_input_buffer(0), std::invalid_argument);
}

TEST(async_input_buffer, buffer_size_huge)
{
    seqan3::dna4_vector vec{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};

    auto v = vec | seqan3::views::async_input_buffer(100000);

    EXPECT_RANGE_EQ(vec, v);
}

TEST(async_input_buffer, destruct_with_full_buffer)
{
    seqan3::dna4_vector vec{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};

    auto v0 = vec | seqan3::views::single_pass_input;

    {
        auto v1 = v0 | seqan3::views::async_input_buffer(5);

        // consume five elements (construction already consumes one)
        auto b = std::ranges::begin(v1);
        ++b; ++b; ++b; ++b;

        /* Give time to rebuffer next five elements so the queue will not be empty.
         * This is not required for this test to be successful, but it is the only
         * way destruction with non-empty buffer is at least likely to happen.
         * And we want it to happen to make sure we don't dead-lock on it.
         */
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    } // thread sync at destruction of v1; tests working destruction with full buffer

    EXPECT_GE(std::ranges::distance(v0), 17); // total of at most 10 chars consumed
}

TEST(async_input_buffer, combinability)
{
    seqan3::dna4_vector vec{"ACGTACGTACGTATCGAGAGCTTTAGC"_dna4};
    seqan3::dna4_vector cmp{"ACGTACGTAC"_dna4};

    auto adapt = seqan3::views::async_input_buffer(5) | seqan3::views::take(10);

    auto v = vec | adapt;

    EXPECT_RANGE_EQ(cmp, v);
}

TEST(async_input_buffer, concepts)
{
    std::vector<int> vec;

    auto v1 = vec | seqan3::views::async_input_buffer(1);

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
}
