// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/concepts>
#include <forward_list>
#include <seqan3/std/iterator>
#include <list>
#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/core/type_traits/iterator.hpp>

template <typename iterator_t>
SEQAN3_CONCEPT iterator_traits_has_iterator_category = requires()
{
    typename std::iterator_traits<iterator_t>::iterator_category;
};

template <typename iterator_t>
SEQAN3_CONCEPT has_iterator_category_tag_t = requires()
{
    typename seqan3::detail::iterator_category_tag_t<iterator_t>;
};

#if SEQAN3_WORKAROUND_GCC_96070
template <typename base_t>
struct my_iterator : base_t
{
    using difference_type = std::iter_difference_t<base_t>;
    using value_type = std::iter_value_t<base_t>;
    using reference = std::iter_reference_t<base_t>;
    using pointer = void;

    using iterator_category = seqan3::detail::iterator_category_tag_t<base_t>;
    using iterator_concept = seqan3::detail::iterator_concept_tag_t<base_t>;
};
#else // ^^^ workaround / no workaround vvv
// See https://github.com/seqan/product_backlog/issues/151 how to make this right, depending on how the standard will
// resolve this issue.
template <typename base_t>
struct inherit_iterator_tag : public base_t
{
    using iterator_concept = seqan3::detail::iterator_concept_tag_t<base_t>;
};

template <has_iterator_category_tag_t base_t>
struct inherit_iterator_tag<base_t> : public base_t
{
    using iterator_category = seqan3::detail::iterator_category_tag_t<base_t>;
    using iterator_concept = seqan3::detail::iterator_concept_tag_t<base_t>;
};

template <typename base_t>
struct my_iterator : inherit_iterator_tag<base_t>
{
    using difference_type = std::iter_difference_t<base_t>;
    using value_type = std::iter_value_t<base_t>;
    using reference = std::iter_reference_t<base_t>;
    using pointer = void;
};
#endif // SEQAN3_WORKAROUND_GCC_96070

#ifdef __cpp_lib_ranges // This is C++20 behaviour.
TEST(iterator_category_tag_t, no_legacy_iterator)
{
    {
        using view_t = std::ranges::basic_istream_view<char, char, std::char_traits<char>>;
        using iterator_t = std::ranges::iterator_t<view_t>;
#if SEQAN3_WORKAROUND_GCC_96070
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>, void>));
#else // ^^^ workaround / no workaround vvv
        EXPECT_FALSE(has_iterator_category_tag_t<iterator_t>);
#endif // SEQAN3_WORKAROUND_GCC_96070
        EXPECT_FALSE(iterator_traits_has_iterator_category<iterator_t>);
    }

    {
        using view_t = std::ranges::basic_istream_view<char, char, std::char_traits<char>>;
        using iterator_t = my_iterator<std::ranges::iterator_t<view_t>>;
#if SEQAN3_WORKAROUND_GCC_96070
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>, void>));
        EXPECT_TRUE((std::same_as<std::iterator_traits<iterator_t>::iterator_category, void>));
#else // ^^^ workaround / no workaround vvv
        EXPECT_FALSE(has_iterator_category_tag_t<iterator_t>);
        EXPECT_FALSE(iterator_traits_has_iterator_category<iterator_t>);
#endif // SEQAN3_WORKAROUND_GCC_96070
    }
}
#endif // __cpp_lib_ranges

TEST(iterator_category_tag_t, output_iterator_tag)
{
    using iterator_t = std::cpp20::ostream_iterator<int>;
    EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                              std::output_iterator_tag>));
}

TEST(iterator_category_tag_t, input_iterator_tag)
{
    {
        using iterator_t = std::istream_iterator<int>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::input_iterator_tag>));
        EXPECT_TRUE((std::same_as<typename my_iterator<iterator_t>::iterator_category,
                                  std::input_iterator_tag>));
    }

#ifdef __cpp_lib_ranges
    {
        // std::views::transform will drop the iterator_category if the lambda doesn't return a lvalue.
        // This is C++20 behaviour.
        using range_t = std::vector<int>;
        auto lambda = [] (auto & element) -> auto { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::input_iterator_tag>));
    }
#endif // __cpp_lib_ranges
}

TEST(iterator_category_tag_t, forward_iterator_tag)
{
    {
        using range_t = std::forward_list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::forward_iterator_tag>));
    }

    {
        // std::views::transform will only keep the iterator_category if the lambda returns a lvalue.
        using range_t = std::forward_list<int>;
        auto lambda = [] (auto & element) -> auto & { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::forward_iterator_tag>));
    }
}

TEST(iterator_category_tag_t, bidirectional_iterator_tag)
{
    {
        using range_t = std::list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::bidirectional_iterator_tag>));
    }

    {
        // std::views::transform will only keep the iterator_category if the lambda returns a lvalue.
        using range_t = std::list<int>;
        auto lambda = [] (auto & element) -> auto & { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::bidirectional_iterator_tag>));
    }
}

TEST(iterator_category_tag_t, random_access_iterator_tag)
{
    {
        using range_t = std::vector<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::random_access_iterator_tag>));
    }

    {
        // std::views::transform will only keep the iterator_category if the lambda returns a lvalue.
        using range_t = std::vector<int>;
        auto lambda = [] (auto & element) -> auto & { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_category_tag_t<iterator_t>,
                                  std::random_access_iterator_tag>));
    }
}

TEST(iterator_concept_tag_t, output_iterator_tag)
{
    using iterator_t = std::cpp20::ostream_iterator<int>;
    EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                              std::output_iterator_tag>));
}

TEST(iterator_concept_tag_t, input_iterator_tag)
{
    {
        using iterator_t = std::istream_iterator<int>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::input_iterator_tag>));
    }

    {
        using view_t = std::ranges::basic_istream_view<char, char, std::char_traits<char>>;
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::input_iterator_tag>));
    }
}

TEST(iterator_concept_tag_t, forward_iterator_tag)
{
    {
        using range_t = std::forward_list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::forward_iterator_tag>));
    }

    {
        // The iterator_concept stays the same, if the lambda returns a lvalue.
        using range_t = std::forward_list<int>;
        auto lambda = [] (auto & element) -> auto & { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::forward_iterator_tag>));
    }

    {
        // Furthermore, the iterator_concept stays the same, even if the lambda doesn't return a lvalue.
        using range_t = std::forward_list<int>;
        auto lambda = [] (auto & element) -> auto { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::forward_iterator_tag>));
    }
}

TEST(iterator_concept_tag_t, bidirectional_iterator_tag)
{
    {
        using range_t = std::list<int>;
        using iterator_t = std::ranges::iterator_t<range_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::bidirectional_iterator_tag>));
    }

    {
        // The iterator_concept stays the same, if the lambda returns a lvalue.
        using range_t = std::list<int>;
        auto lambda = [] (auto & element) -> auto & { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::bidirectional_iterator_tag>));
    }

    {
        // Furthermore, the iterator_concept stays the same, even if the lambda doesn't return a lvalue.
        using range_t = std::list<int>;
        auto lambda = [] (auto & element) -> auto { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::bidirectional_iterator_tag>));
    }
}

TEST(iterator_concept_tag_t, random_access_iterator_tag)
{
    {
        // A contiguous_iterator will be downgraded to a random_access_iterator by std::views::transform
        using range_t = std::vector<int>; // contiguous_iterator
        auto lambda = [] (auto & element) -> auto & { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::random_access_iterator_tag>));
    }

    {
        // Furthermore, the iterator_concept stays the same, even if the lambda doesn't return a lvalue.
        using range_t = std::vector<int>; // contiguous_iterator
        auto lambda = [] (auto & element) -> auto { return element; };
        using view_t = decltype(std::declval<range_t &>() | std::views::transform(lambda));
        using iterator_t = std::ranges::iterator_t<view_t>;
        EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                                  std::random_access_iterator_tag>));
    }
}

TEST(iterator_concept_tag_t, contiguous_iterator_tag)
{
    using range_t = std::vector<int>;
    using iterator_t = std::ranges::iterator_t<range_t>;
    EXPECT_TRUE((std::same_as<seqan3::detail::iterator_concept_tag_t<iterator_t>,
                              std::contiguous_iterator_tag>));
}
