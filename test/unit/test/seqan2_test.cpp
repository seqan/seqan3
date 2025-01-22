// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ranges>

#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/test/seqan2.hpp>

// Note: this file will only test regressions encountered with seqan2 compatibility and has no claim to be complete

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/sequence.h>

template <typename T>
class seqan2_container : public ::testing::Test
{};

using seqan2_container_types = ::testing::Types<seqan2::String<int>, seqan2::StringSet<int>>;
TYPED_TEST_SUITE(seqan2_container, seqan2_container_types, );

template <typename container_t>
container_t construct_iota(int n)
{
    container_t container{};

    for (int i = 0; i < n; ++i)
        seqan2::appendValue(container, i);

    return container;
}

TYPED_TEST(seqan2_container, std_ranges_size)
{
    using container_t = TypeParam;
    container_t container = construct_iota<container_t>(5u);

    EXPECT_EQ(5u, std::ranges::size(container));
}

TYPED_TEST(seqan2_container, std_ranges_begin_end)
{
    using container_t = TypeParam;
    container_t container = construct_iota<container_t>(5u);

    auto it = std::ranges::begin(container);
    auto it_end = std::ranges::end(container);

    EXPECT_SAME_TYPE(decltype(it), decltype(seqan2::begin(container)));
    EXPECT_SAME_TYPE(decltype(it_end), decltype(seqan2::end(container)));

    for (int i = 0; it != it_end; ++it, ++i)
    {
        EXPECT_EQ(i, *it);
    }
}

TYPED_TEST(seqan2_container, std_ranges_iterator)
{
    using container_t = TypeParam;
    using iterator_t = decltype(std::ranges::begin(std::declval<container_t &>()));
    using const_iterator_t = decltype(std::ranges::begin(std::declval<container_t const &>()));

    EXPECT_SAME_TYPE(std::ranges::iterator_t<container_t>, iterator_t);
    EXPECT_SAME_TYPE(std::ranges::iterator_t<container_t &>, iterator_t);
    EXPECT_SAME_TYPE(std::ranges::iterator_t<container_t const>, const_iterator_t);
    EXPECT_SAME_TYPE(std::ranges::iterator_t<container_t const &>, const_iterator_t);
}

TYPED_TEST(seqan2_container, std_iterator_traits)
{
    using container_t = TypeParam;
    using iterator_t = std::ranges::iterator_t<container_t>;
    using value_type = typename std::iterator_traits<iterator_t>::value_type;
    EXPECT_SAME_TYPE(value_type, int);
}

TYPED_TEST(seqan2_container, std_iterator)
{
    using container_t = TypeParam;
    using iterator_t = std::ranges::iterator_t<container_t>;
    EXPECT_FALSE(std::input_or_output_iterator<container_t>);
    EXPECT_TRUE(std::input_or_output_iterator<iterator_t>);
}

template <typename range_t>
concept SeqAn2Range = requires (range_t range) {
    { seqan2::begin(range) } -> std::same_as<std::ranges::iterator_t<range_t>>;
    { seqan2::end(range) } -> std::same_as<std::ranges::iterator_t<range_t>>;
};

TYPED_TEST(seqan2_container, seqan_range_concept)
{
    using container_t = TypeParam;
    using iterator_t = std::ranges::iterator_t<container_t>;
    EXPECT_TRUE(SeqAn2Range<container_t>);
    EXPECT_FALSE(SeqAn2Range<iterator_t>);
}

TYPED_TEST(seqan2_container, std_ranges_range)
{
    using container_t = TypeParam;
    using iterator_t = std::ranges::iterator_t<container_t>;
    EXPECT_TRUE(std::ranges::range<container_t>);
    EXPECT_FALSE(std::ranges::range<iterator_t>);
}

TYPED_TEST(seqan2_container, seqan3_value_type)
{
    using container_t = TypeParam;
    using value_type = std::ranges::range_value_t<container_t>;
    EXPECT_SAME_TYPE(value_type, int);
}

#endif
