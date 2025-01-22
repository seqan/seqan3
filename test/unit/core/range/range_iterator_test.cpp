// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <list>
#include <ranges>
#include <vector>

#include <seqan3/core/range/detail/random_access_iterator.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/type_list/traits.hpp>

TEST(range_and_iterator, iterator_t)
{
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<std::vector<int>>, typename std::vector<int>::iterator>));
    EXPECT_TRUE(
        (std::is_same_v<std::ranges::iterator_t<std::vector<int> const>, typename std::vector<int>::const_iterator>));

    auto v = std::views::iota(1);
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<decltype(v)>, decltype(std::ranges::begin(v))>));
    EXPECT_FALSE((std::is_same_v<std::ranges::iterator_t<decltype(v)>, decltype(std::ranges::end(v))>));
}

TEST(range_and_iterator, maybe_const_iterator_t)
{
    // common range
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<false, std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, std::vector<int>>,
                                typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<false, std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));

    // non-common range
    auto v = std::views::iota(1);
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::maybe_const_iterator_t<false, decltype(v)>, decltype(std::ranges::begin(v))>));
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, decltype(v)>, decltype(std::ranges::cbegin(v))>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<false, decltype(v) const>,
                                decltype(std::ranges::cbegin(v))>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, decltype(v) const>,
                                decltype(std::ranges::cbegin(v))>));

    EXPECT_FALSE(
        (std::is_same_v<seqan3::detail::maybe_const_iterator_t<false, decltype(v)>, decltype(std::ranges::end(v))>));
    EXPECT_FALSE(
        (std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, decltype(v)>, decltype(std::ranges::cend(v))>));
    EXPECT_FALSE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, decltype(v) const>,
                                 decltype(std::ranges::cend(v))>));
    EXPECT_FALSE((std::is_same_v<seqan3::detail::maybe_const_iterator_t<true, decltype(v) const>,
                                 decltype(std::ranges::cend(v))>));
}

TEST(range_and_iterator, sentinel_t)
{
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<std::vector<int>>, typename std::vector<int>::iterator>));
    EXPECT_TRUE(
        (std::is_same_v<std::ranges::sentinel_t<std::vector<int> const>, typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<std::vector<int>>, std::ranges::iterator_t<std::vector<int>>>));

    auto v = std::views::iota(1);
    EXPECT_FALSE((std::is_same_v<std::ranges::sentinel_t<decltype(v)>, decltype(std::ranges::begin(v))>));
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<decltype(v)>, decltype(std::ranges::end(v))>));
}

TEST(range_and_iterator, maybe_const_sentinel_t)
{
    // common range
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<false, std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<true, std::vector<int>>,
                                typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<false, std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<true, std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));

    // non-common range
    auto v = std::views::iota(1);
    EXPECT_FALSE(
        (std::is_same_v<seqan3::detail::maybe_const_sentinel_t<false, decltype(v)>, decltype(std::ranges::begin(v))>));
    EXPECT_FALSE(
        (std::is_same_v<seqan3::detail::maybe_const_sentinel_t<true, decltype(v)>, decltype(std::ranges::cbegin(v))>));
    EXPECT_FALSE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<false, decltype(v) const>,
                                 decltype(std::ranges::cbegin(v))>));
    EXPECT_FALSE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<true, decltype(v) const>,
                                 decltype(std::ranges::cbegin(v))>));

    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::maybe_const_sentinel_t<false, decltype(v)>, decltype(std::ranges::end(v))>));
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::maybe_const_sentinel_t<true, decltype(v)>, decltype(std::ranges::cend(v))>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<false, decltype(v) const>,
                                decltype(std::ranges::cend(v))>));
    EXPECT_TRUE((std::is_same_v<seqan3::detail::maybe_const_sentinel_t<true, decltype(v) const>,
                                decltype(std::ranges::cend(v))>));
}

template <typename list1, typename list2, size_t pos = 0>
void expect_same_types()
{
    constexpr bool val = std::is_same_v<seqan3::list_traits::at<pos, list1>, seqan3::list_traits::at<pos, list2>>;
    if constexpr (!val)
    {
        std::cerr << "pos: " << pos << " \'" << seqan3::detail::type_name_as_string<seqan3::list_traits::at<pos, list1>>
                  << "\' not equal to \'" << seqan3::detail::type_name_as_string<seqan3::list_traits::at<pos, list2>>
                  << "\' \n";
    }
    EXPECT_TRUE(val);

    if constexpr (pos < list1::size() - 1)
        expect_same_types<list1, list2, pos + 1>();
}

TEST(range_and_iterator, range_iter_value_t)
{
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);

    using type_list_example = seqan3::type_list<std::ranges::range_value_t<std::vector<int>>,       // short
                                                typename std::vector<int>::value_type,              // member type
                                                std::ranges::range_value_t<std::vector<int> const>, // const container
                                                std::iter_value_t<iterator_of_int_vector>,          // iterator
                                                std::iter_value_t<foreign_iterator>,                // iterator2
                                                std::ranges::range_value_t<decltype(v)>>;           // range, no member

    using comp_list = seqan3::type_list<int, int, int, int, int, int>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, range_iter_reference_t)
{
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example =
        seqan3::type_list<std::ranges::range_reference_t<std::vector<int>>,       // short
                          typename std::vector<int>::reference,                   // member type
                          std::ranges::range_reference_t<std::vector<int> const>, // const container
                          std::iter_reference_t<iterator_of_int_vector>,          // iterator
                          std::iter_reference_t<foreign_iterator>,                // iterator2
                          std::ranges::range_reference_t<decltype(v)>>;           // range, no member

    using comp_list = seqan3::type_list<int &,
                                        int &,
                                        int const &, // container is const
                                        int &,
                                        int &,
                                        int>; // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, range_iter_rvalue_reference_t)
{
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = seqan3::type_list<
        std::ranges::range_rvalue_reference_t<std::vector<int>>, // short
                                                                 // No types have member,yet:
        //                                              typename std::vector<int>::rvalue_reference, // member type
        std::ranges::range_rvalue_reference_t<std::vector<int> const>, // const container
        std::iter_rvalue_reference_t<iterator_of_int_vector>,          // iterator
        std::iter_rvalue_reference_t<foreign_iterator>,                // iterator2
        std::ranges::range_rvalue_reference_t<decltype(v)>>;           // range, no member

    using comp_list = seqan3::type_list<int &&,
                                        //                                      int &&,
                                        int const &&, // container is const
                                        int &&,
                                        int &&,
                                        int>; // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, const_range_iter_reference_t)
{
    using const_iterator_of_int_vector = std::ranges::iterator_t<std::vector<int> const>;
    using const_foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int> const>;
    using const_iterator_of_view = std::ranges::iterator_t<decltype(std::views::iota(1)) const>;
    using type_list_example = seqan3::type_list<std::ranges::range_reference_t<std::vector<int> const>, // short
                                                typename std::vector<int>::const_reference,             // member type
                                                std::iter_reference_t<const_iterator_of_int_vector>,    // iterator
                                                std::iter_reference_t<const_foreign_iterator>,          // iterator2
                                                std::iter_reference_t<const_iterator_of_view>>; // range, no member

    using comp_list = seqan3::type_list<int const &, int const &, int const &, int const &,
                                        int>; // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, range_iter_difference_t)
{
    //TODO(h-2): add something that actually has a different difference_type
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    using type_list_example =
        seqan3::type_list<std::ranges::range_difference_t<std::vector<int>>,       // short
                          typename std::vector<int>::difference_type,              // member type
                          std::ranges::range_difference_t<std::vector<int> const>, // const container
                          std::iter_difference_t<iterator_of_int_vector>,          // iterator
                          std::iter_difference_t<foreign_iterator>>;               // range, no member

    using comp_list = seqan3::type_list<std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t, std::ptrdiff_t>;

    expect_same_types<type_list_example, comp_list>();

    // views::ints' difference_type is not std::ptrdiff_t, but depends on the size.
    // For infinite views, like in our case, it's std::int_fast64_t (or std::int_fast32_t on 32bit).
    // On most platforms this is the same as std::size_t (`long int` on Linux and FreeBSD;
    // `long long int` on macOS and MSVC) BUT on some platforms like OpenBSD std::ptrdiff_t is
    // `long int`, but std::int64_t and std::int_fast64_t are `long long int`.
    // This detects word-size and selects the corresponding int_fast*_t which is always correct
    using view_int_diff_t = std::conditional_t<sizeof(std::ptrdiff_t) == 8, std::int_fast64_t, std::int_fast32_t>;

    auto v = std::views::iota(1);
    using iota_difference_t = std::ranges::range_difference_t<decltype(v)>;
    EXPECT_TRUE((std::signed_integral<iota_difference_t>));
    EXPECT_EQ(sizeof(iota_difference_t), sizeof(view_int_diff_t));
}

TEST(range_and_iterator, range_size_t)
{
    //TODO(h-2): add something that actually has a different size_type
    // iota is not a sized range, but take_exactly is
    auto v = std::views::iota(0) | seqan3::detail::take_exactly(2);
    using type_list_example = seqan3::type_list<std::ranges::range_size_t<std::vector<int>>,       // short
                                                typename std::vector<int>::size_type,              // member type
                                                std::ranges::range_size_t<std::vector<int> const>, // const container
                                                std::ranges::range_size_t<decltype(v)>>;           // range, no member

    using comp_list = seqan3::type_list<size_t, size_t, size_t, size_t>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, range_innermost_value)
{
    using vector_of_int_vector = std::vector<std::vector<int>>;
    using type_list_example = seqan3::type_list<typename seqan3::range_innermost_value<std::vector<int>>::type, // long
                                                seqan3::range_innermost_value_t<std::vector<int>>,              // short
                                                seqan3::range_innermost_value_t<vector_of_int_vector>>; // two-level

    using comp_list = seqan3::type_list<int, int, int, int, int>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, range_dimension_v)
{
    EXPECT_EQ(1u, seqan3::range_dimension_v<std::vector<int>>);
    EXPECT_EQ(2u, seqan3::range_dimension_v<std::vector<std::vector<int>>>);
}
