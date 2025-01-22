// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>
#include <list>
#include <ranges>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/repeat.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

#include "../../range/iterator_test_template.hpp"

using iterator_type = std::ranges::iterator_t<decltype(std::declval<std::vector<int> &>() | seqan3::views::chunk(4))>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<int> text{1, 4, 2, 7, 4, 5, 8, 3, 4, 7, 5, 4, 3};
    std::vector<std::vector<int>> expected_range{{1, 4, 2, 7}, {4, 5, 8, 3}, {4, 7, 5, 4}, {3}};

    using test_range_t = decltype(text | seqan3::views::chunk(4));
    test_range_t test_range = text | seqan3::views::chunk(4);

    template <typename lhs_t, typename rhs_t>
    static void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        EXPECT_RANGE_EQ(lhs, rhs);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, iterator_type, );

using underlying_range_types = ::testing::Types<std::forward_list<int>,       // forward range
                                                std::forward_list<int> const, // forward range
                                                std::list<int>,               // bidirectional range
                                                std::list<int> const,         // bidirectional range
                                                std::vector<int>,             // random access range
                                                std::vector<int> const>;      // random access range

template <class range_type>
struct chunk_view_test : public ::testing::Test
{
    range_type text{1, 4, 10, 2, 7};
    decltype(text | seqan3::views::chunk(2)) v = text | seqan3::views::chunk(2);
};

TYPED_TEST_SUITE(chunk_view_test, underlying_range_types, );

TYPED_TEST(chunk_view_test, concepts)
{
    // chunk view preserves input / forward / bidirectional / random_access / sized / common.
    EXPECT_EQ((std::ranges::range<TypeParam>), (std::ranges::range<decltype(this->v)>));
    EXPECT_EQ((std::ranges::input_range<TypeParam>), (std::ranges::input_range<decltype(this->v)>));
    EXPECT_EQ((std::ranges::forward_range<TypeParam>), (std::ranges::forward_range<decltype(this->v)>));
    EXPECT_EQ((std::ranges::bidirectional_range<TypeParam>), (std::ranges::bidirectional_range<decltype(this->v)>));
    EXPECT_EQ((std::ranges::random_access_range<TypeParam>), (std::ranges::random_access_range<decltype(this->v)>));
    EXPECT_EQ((std::ranges::sized_range<TypeParam>), (std::ranges::sized_range<decltype(this->v)>));
    EXPECT_EQ((std::ranges::common_range<TypeParam>), (std::ranges::common_range<decltype(this->v)>));

    // it always ensures view
    EXPECT_TRUE((std::ranges::view<decltype(this->v)>));

    // it loses contiguous range and output range
    EXPECT_FALSE((std::ranges::contiguous_range<decltype(this->v)>));
    EXPECT_FALSE((std::ranges::output_range<decltype(this->v), std::ranges::range_value_t<TypeParam> &>));
}

TYPED_TEST(chunk_view_test, construction)
{
    EXPECT_EQ(std::default_initializable<decltype(this->v)>,
              std::default_initializable<std::views::all_t<decltype(this->text) &>>);
    EXPECT_TRUE((std::is_copy_constructible_v<decltype(this->v)>));
    EXPECT_TRUE((std::is_move_constructible_v<decltype(this->v)>));
    EXPECT_TRUE((std::is_copy_assignable_v<decltype(this->v)>));
    EXPECT_TRUE((std::is_move_assignable_v<decltype(this->v)>));
}

TYPED_TEST(chunk_view_test, distance_and_size)
{
    EXPECT_EQ(std::ranges::distance(this->v), 3);

    if constexpr (std::ranges::sized_range<TypeParam>)
    {
        EXPECT_EQ(this->v.size(), 3u);
    }
}

TYPED_TEST(chunk_view_test, view_compatibility_test)
{
    {
        auto v = this->text
               | std::views::transform(
                     [](int i)
                     {
                         return i + 1;
                     })
               | seqan3::views::chunk(2);

        std::vector<std::vector<int>> expected{{2, 5}, {11, 3}, {8}};

        auto v_it = v.begin();
        auto expected_it = expected.begin();
        for (; v_it != v.end(); ++v_it, ++expected_it)
            EXPECT_RANGE_EQ(*v_it, *expected_it);
    }

    {
        auto v = this->text | seqan3::views::chunk(2)
               | std::views::transform(
                     [](auto chunk)
                     {
                         return std::ranges::distance(chunk);
                     });

        std::vector<size_t> expected{2, 2, 1};

        EXPECT_RANGE_EQ(v, expected);
    }

    { // combine with a range that does not model std::ranges::common_range
        ASSERT_FALSE((std::ranges::common_range<decltype(seqan3::views::repeat(42))>)); // sanity check

        auto v = seqan3::views::repeat(42) | seqan3::views::chunk(2);

        std::vector<std::vector<int>> expected{{42, 42}, {42, 42}};

        auto v_it = v.begin();
        auto expected_it = expected.begin();
        for (size_t i = 0; i < 2; ++i, ++v_it, ++expected_it)
            EXPECT_RANGE_EQ(*v_it, *expected_it);
    }
}

TYPED_TEST(chunk_view_test, underlying_input_range_test)
{
    { // fully consume chunks directly
        auto v = this->text | seqan3::views::single_pass_input | seqan3::views::chunk(2);

        std::vector<std::vector<int>> expected{{1, 4}, {10, 2}, {7}};

        auto v_it = v.begin();
        auto expected_it = expected.begin();
        for (; v_it != v.end(); ++v_it, ++expected_it)
            EXPECT_RANGE_EQ(*v_it, *expected_it);
    }

    { // do not fully consume chunks directly but let the iterator consume the chunk
        auto v = this->text | seqan3::views::single_pass_input | seqan3::views::chunk(2);

        std::vector<int> expected{1, 10, 7};

        auto v_it = v.begin();
        auto expected_it = expected.begin();
        for (; v_it != v.end(); ++v_it, ++expected_it)
            EXPECT_EQ(*(*v_it).begin(), *expected_it);
    }
}

TYPED_TEST(chunk_view_test, use_on_temporaries)
{
    if constexpr (!std::is_const_v<TypeParam>)
    {
        std::vector<std::vector<int>> expected_range{{1, 4, 2, 7}, {4, 5, 8, 3}, {4, 7, 5, 4}, {3}};

        size_t i{};
        for (auto && chunk : seqan3::views::chunk(TypeParam{1, 4, 2, 7, 4, 5, 8, 3, 4, 7, 5, 4, 3}, 4))
        {
            EXPECT_RANGE_EQ(chunk, expected_range[i]);
            ++i;
        }
        EXPECT_EQ(i, 4u);
    }
}
