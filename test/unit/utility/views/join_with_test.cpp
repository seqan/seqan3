// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>
#include <list>
#include <vector>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/join_with.hpp>
#include <seqan3/utility/views/repeat_n.hpp>

#include "../../range/iterator_test_template.hpp"

using range_t = std::vector<std::string>;
using pattern_t = std::vector<char>;
using join_with_view_t = decltype(seqan3::views::join_with(std::declval<range_t &>(), std::declval<pattern_t &>()));
using join_with_iterator_t = std::ranges::iterator_t<join_with_view_t>;

template <>
struct iterator_fixture<join_with_iterator_t> : public ::testing::Test
{
    using iterator_tag = std::bidirectional_iterator_tag;

    static constexpr bool const_iterable = true;

    range_t range{"AA", "BBB", "CC", "DDD"};

    pattern_t pattern{'x', 'y'};

    std::string expected_range{"AAxyBBBxyCCxyDDD"};

    join_with_view_t test_range{seqan3::views::join_with(range, pattern)};
};
INSTANTIATE_TYPED_TEST_SUITE_P(join_with_iterator_test, iterator_fixture, join_with_iterator_t, );

class join_with : public ::testing::Test
{
protected:
    using forward_range_t = std::forward_list<std::vector<int>>;
    using bidirectional_range_t = std::list<std::vector<int>>;
    using random_access_range_t = std::vector<std::vector<int>>;

    using forward_pattern_t = std::forward_list<int>;
    using bidirectional_pattern_t = std::list<int>;
    using random_access_pattern_t = std::vector<int>;

    forward_range_t forward_range{{0, 1}, {2, 3}, {3, 4}};
    bidirectional_range_t bidirectional_range{forward_range.begin(), forward_range.end()};
    random_access_range_t random_access_range{forward_range.begin(), forward_range.end()};

    forward_pattern_t forward_pattern{23};
    bidirectional_pattern_t bidirectional_pattern{forward_pattern.begin(), forward_pattern.end()};
    random_access_pattern_t random_access_pattern{forward_pattern.begin(), forward_pattern.end()};

    using join_with_forward_1_t = decltype(seqan3::views::join_with(forward_range, forward_pattern));
    using join_with_forward_2_t = decltype(seqan3::views::join_with(forward_range, bidirectional_pattern));
    using join_with_forward_3_t = decltype(seqan3::views::join_with(bidirectional_range, forward_pattern));

    using join_with_bidirectional_1_t = decltype(seqan3::views::join_with(bidirectional_range, bidirectional_pattern));
    using join_with_bidirectional_2_t = decltype(seqan3::views::join_with(bidirectional_range, random_access_pattern));
    using join_with_bidirectional_3_t = decltype(seqan3::views::join_with(random_access_range, bidirectional_pattern));

    using join_with_random_access_t = decltype(seqan3::views::join_with(random_access_range, random_access_pattern));
};

TEST_F(join_with, concepts)
{
    // {forward, bidirectional, random_access}_range
    EXPECT_TRUE(std::ranges::forward_range<join_with_forward_1_t>);
    EXPECT_FALSE(std::ranges::bidirectional_range<join_with_forward_1_t>);
    EXPECT_TRUE(std::ranges::forward_range<join_with_forward_2_t>);
    EXPECT_FALSE(std::ranges::bidirectional_range<join_with_forward_2_t>);
    EXPECT_TRUE(std::ranges::forward_range<join_with_forward_3_t>);
    EXPECT_FALSE(std::ranges::bidirectional_range<join_with_forward_3_t>);

    EXPECT_TRUE(std::ranges::bidirectional_range<join_with_bidirectional_1_t>);
    EXPECT_FALSE(std::ranges::random_access_range<join_with_bidirectional_1_t>);
    EXPECT_TRUE(std::ranges::bidirectional_range<join_with_bidirectional_2_t>);
    EXPECT_FALSE(std::ranges::random_access_range<join_with_bidirectional_2_t>);
    EXPECT_TRUE(std::ranges::bidirectional_range<join_with_bidirectional_3_t>);
    EXPECT_FALSE(std::ranges::random_access_range<join_with_bidirectional_3_t>);

    EXPECT_FALSE(std::ranges::random_access_range<join_with_random_access_t>);

    // common_range
    EXPECT_TRUE(std::ranges::common_range<join_with_forward_1_t>);
    EXPECT_TRUE(std::ranges::common_range<join_with_forward_2_t>);
    EXPECT_TRUE(std::ranges::common_range<join_with_forward_3_t>);

    EXPECT_TRUE(std::ranges::common_range<join_with_bidirectional_1_t>);
    EXPECT_TRUE(std::ranges::common_range<join_with_bidirectional_2_t>);
    EXPECT_TRUE(std::ranges::common_range<join_with_bidirectional_3_t>);

    EXPECT_TRUE(std::ranges::common_range<join_with_random_access_t>);

    // // No common_range because the range_reference_t must be a reference. repeat_n returns values.
    auto repeat_n_view = seqan3::views::repeat_n(2, 2);
    std::vector<decltype(repeat_n_view)> test_range{repeat_n_view, repeat_n_view};
    EXPECT_FALSE(std::ranges::common_range<decltype(seqan3::views::join_with(test_range, forward_pattern))>);
}

TEST_F(join_with, basic)
{
    {
        std::vector<int> const expected_result{0, 1, 23, 2, 3, 23, 3, 4};
        EXPECT_RANGE_EQ(seqan3::views::join_with(forward_range, forward_pattern), expected_result);
    }
    {
        std::vector<int> const expected_result{2, 2, 23, 2, 2};
        auto repeat_n_view = seqan3::views::repeat_n(2, 2);
        std::vector<decltype(repeat_n_view)> test_range{repeat_n_view, repeat_n_view};
        EXPECT_RANGE_EQ(seqan3::views::join_with(test_range, forward_pattern), expected_result);
    }
}

TEST_F(join_with, combine)
{
    {
        std::vector<int> const expected_result{4, 3, 23, 3, 2, 23, 1, 0};
        EXPECT_RANGE_EQ(seqan3::views::join_with(bidirectional_range, bidirectional_pattern) | std::views::reverse,
                        expected_result);
    }
    {
        std::vector<int> const expected_result{2, 2, 23, 2, 2};
        auto repeat_n_view = seqan3::views::repeat_n(2, 2);
        std::vector<decltype(repeat_n_view)> test_range{repeat_n_view, repeat_n_view};
        EXPECT_RANGE_EQ(seqan3::views::join_with(test_range, bidirectional_pattern) | std::views::all, expected_result);
    }
}
