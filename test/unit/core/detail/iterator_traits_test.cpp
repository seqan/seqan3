// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <concepts>
#include <forward_list>
#include <iterator>
#include <list>
#include <ranges>
#include <vector>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/test/expect_same_type.hpp>

template <typename t>
concept has_iterator_category = requires () { typename t::iterator_category; };

template <typename it_t>
struct iterator_category_tag_or_void
{
    using type = void;
};

template <typename it_t>
    requires has_iterator_category<std::iterator_traits<it_t>>
struct iterator_category_tag_or_void<it_t>
{
    using type = typename std::iterator_traits<it_t>::iterator_category;
};

template <typename it_t>
using iterator_category_tag_or_void_t = typename iterator_category_tag_or_void<it_t>::type;

template <typename base_t>
struct my_iterator : public base_t, public seqan3::detail::maybe_inherited_iterator_category<base_t>
{
    using difference_type = std::iter_difference_t<base_t>;
    using value_type = std::iter_value_t<base_t>;
    using reference = std::iter_reference_t<base_t>;
    using pointer = void;
    using iterator_concept = seqan3::detail::iterator_concept_tag_t<base_t>;
};

TEST(maybe_iterator_category, no_legacy_iterator)
{
    {
        using view_t = std::ranges::basic_istream_view<char, char, std::char_traits<char>>;
        using iterator_t = std::ranges::iterator_t<view_t>;
        // not defined (this is expected for C++20 input iterator)
        // EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, void);
        EXPECT_FALSE(has_iterator_category<iterator_t>);
        EXPECT_FALSE(has_iterator_category<std::iterator_traits<iterator_t>>);
    }

    {
        using view_t = std::ranges::basic_istream_view<char, char, std::char_traits<char>>;
        using iterator_t = my_iterator<std::ranges::iterator_t<view_t>>;
        EXPECT_FALSE(has_iterator_category<iterator_t>);
        EXPECT_FALSE(has_iterator_category<std::iterator_traits<iterator_t>>);
    }
}

TEST(maybe_iterator_category, output_iterator_tag)
{
    {
        using iterator_t = std::ostream_iterator<int>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::output_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::output_iterator_tag);
    }

    {
        using iterator_t = std::ostream_iterator<int>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::output_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::output_iterator_tag);
    }
}

TEST(maybe_iterator_category, input_iterator_tag)
{
    {
        using iterator_t = std::istream_iterator<int>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::input_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::input_iterator_tag);
    }

    {
        // std::views::transform will drop the iterator_category if the lambda doesn't return a lvalue.
        // This is C++20 behaviour.
        using range_t = std::vector<int>;
        auto lambda = [](auto & element) -> auto
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::input_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::input_iterator_tag);
    }
}

TEST(maybe_iterator_category, forward_iterator_tag)
{
    {
        using range_t = std::forward_list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::forward_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::forward_iterator_tag);
    }

    {
        // std::views::transform will only keep the iterator_category if the lambda returns a lvalue.
        using range_t = std::forward_list<int>;
        auto lambda = [](auto & element) -> auto &
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::forward_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::forward_iterator_tag);
    }
}

TEST(maybe_iterator_category, bidirectional_iterator_tag)
{
    {
        using range_t = std::list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::bidirectional_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::bidirectional_iterator_tag);
    }

    {
        // std::views::transform will only keep the iterator_category if the lambda returns a lvalue.
        using range_t = std::list<int>;
        auto lambda = [](auto & element) -> auto &
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::bidirectional_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::bidirectional_iterator_tag);
    }
}

TEST(maybe_iterator_category, random_access_iterator_tag)
{
    {
        using range_t = std::vector<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::random_access_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::random_access_iterator_tag);
    }

    {
        // std::views::transform will only keep the iterator_category if the lambda returns a lvalue.
        using range_t = std::vector<int>;
        auto lambda = [](auto & element) -> auto &
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(std::iterator_traits<iterator_t>::iterator_category, std::random_access_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_category, std::random_access_iterator_tag);
    }
}

template <typename t>
concept has_iterator_concept = requires () { typename t::iterator_concept; };

TEST(iterator_concept_tag_t, output_iterator_tag)
{
    using iterator_t = std::ostream_iterator<int>;
    EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::output_iterator_tag);
    EXPECT_FALSE(has_iterator_concept<iterator_t>);
}

TEST(iterator_concept_tag_t, input_iterator_tag)
{
    {
        using iterator_t = std::istream_iterator<int>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::input_iterator_tag);
        EXPECT_FALSE(has_iterator_concept<iterator_t>);
    }

    {
        using view_t = std::ranges::basic_istream_view<char, char, std::char_traits<char>>;
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::input_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::input_iterator_tag);
    }
}

TEST(iterator_concept_tag_t, forward_iterator_tag)
{
    {
        using range_t = std::forward_list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::forward_iterator_tag);
        EXPECT_FALSE(has_iterator_concept<iterator_t>);
    }

    {
        // The iterator_concept stays the same, if the lambda returns a lvalue.
        using range_t = std::forward_list<int>;
        auto lambda = [](auto & element) -> auto &
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::forward_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::forward_iterator_tag);
    }

    {
        // Furthermore, the iterator_concept stays the same, even if the lambda doesn't return a lvalue.
        using range_t = std::forward_list<int>;
        auto lambda = [](auto & element) -> auto
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::forward_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::forward_iterator_tag);
    }
}

TEST(iterator_concept_tag_t, bidirectional_iterator_tag)
{
    {
        using range_t = std::list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::bidirectional_iterator_tag);
        EXPECT_FALSE(has_iterator_concept<iterator_t>);
    }

    {
        // The iterator_concept stays the same, if the lambda returns a lvalue.
        using range_t = std::list<int>;
        auto lambda = [](auto & element) -> auto &
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::bidirectional_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::bidirectional_iterator_tag);
    }

    {
        // Furthermore, the iterator_concept stays the same, even if the lambda doesn't return a lvalue.
        using range_t = std::list<int>;
        auto lambda = [](auto & element) -> auto
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::bidirectional_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::bidirectional_iterator_tag);
    }
}

TEST(iterator_concept_tag_t, random_access_iterator_tag)
{
    {
        // A contiguous_iterator will be downgraded to a random_access_iterator by std::views::transform
        using range_t = std::vector<int>; // contiguous_iterator
        auto lambda = [](auto & element) -> auto &
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::random_access_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::random_access_iterator_tag);
    }

    {
        // Furthermore, the iterator_concept stays the same, even if the lambda doesn't return a lvalue.
        using range_t = std::vector<int>; // contiguous_iterator
        auto lambda = [](auto & element) -> auto
        {
            return element;
        };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::random_access_iterator_tag);
        EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::random_access_iterator_tag);
    }
}

TEST(iterator_concept_tag_t, contiguous_iterator_tag)
{
    using range_t = std::vector<int>;
    using iterator_t = std::ranges::iterator_t<range_t>;
    EXPECT_SAME_TYPE(seqan3::detail::iterator_concept_tag_t<iterator_t>, std::contiguous_iterator_tag);
    EXPECT_SAME_TYPE(iterator_t::iterator_concept, std::contiguous_iterator_tag);
}
