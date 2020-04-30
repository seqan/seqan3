// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <list>
#include <vector>

#include <gtest/gtest.h>

#include <meta/meta.hpp>


#include <seqan3/core/type_list/all.hpp>
#include <seqan3/core/type_traits/all.hpp>
#include <seqan3/core/detail/type_inspection.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/std/ranges>

TEST(range_and_iterator, iterator_)
{
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));

    auto v = std::views::iota(1);
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<decltype(v)>,
                                decltype(std::ranges::begin(v))>));
    EXPECT_FALSE((std::is_same_v<std::ranges::iterator_t<decltype(v)>,
                                 decltype(std::ranges::end(v))>));
}

TEST(range_and_iterator, sentinel_)
{
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<std::vector<int>>,
                                std::ranges::iterator_t<std::vector<int>>>));

    auto v = std::views::iota(1);
    EXPECT_FALSE((std::is_same_v<std::ranges::sentinel_t<decltype(v)>,
                                 decltype(std::ranges::begin(v))>));
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<decltype(v)>,
                                decltype(std::ranges::end(v))>));
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

TEST(range_and_iterator, value_type_)
{
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);

    using type_list_example = seqan3::type_list<std::ranges::range_value_t<std::vector<int>>, // short
                                                typename std::vector<int>::value_type, // member type
                                                std::ranges::range_value_t<std::vector<int> const>, // const container
                                                std::iter_value_t<iterator_of_int_vector>, // iterator
                                                std::iter_value_t<foreign_iterator>, // iterator2
                                                std::ranges::range_value_t<decltype(v)>>; // range, no member

    using comp_list = seqan3::type_list<int,
                                        int,
                                        int,
                                        int,
                                        int,
                                        int>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, reference_)
{
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = seqan3::type_list<std::ranges::range_reference_t<std::vector<int>>, // short
                                                typename std::vector<int>::reference, // member type
                                                std::ranges::range_reference_t<std::vector<int> const>, // const container
                                                std::iter_reference_t<iterator_of_int_vector>, // iterator
                                                std::iter_reference_t<foreign_iterator>, // iterator2
                                                std::ranges::range_reference_t<decltype(v)>>; // range, no member

    using comp_list = seqan3::type_list<int &,
                                        int &,
                                        int const &, // container is const
                                        int &,
                                        int &,
                                        int>; // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, rvalue_reference_)
{
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = seqan3::type_list<std::ranges::range_rvalue_reference_t<std::vector<int>>, // short
// No types have member,yet:
//                                              typename std::vector<int>::rvalue_reference, // member type
                                                std::ranges::range_rvalue_reference_t<std::vector<int> const>, // const container
                                                std::iter_rvalue_reference_t<iterator_of_int_vector>, // iterator
                                                std::iter_rvalue_reference_t<foreign_iterator>, // iterator2
                                                std::ranges::range_rvalue_reference_t<decltype(v)>>; // range, no member

    using comp_list = seqan3::type_list<int &&,
//                                      int &&,
                                        int const &&,  // container is const
                                        int &&,
                                        int &&,
                                        int>; // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, const_reference_)
{
    using const_iterator_of_int_vector = std::ranges::iterator_t<std::vector<int> const>;
    using const_foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int> const>;
    using const_iterator_of_view = std::ranges::iterator_t<decltype(std::views::iota(1)) const>;
    using type_list_example = seqan3::type_list<std::ranges::range_reference_t<std::vector<int> const>, // short
                                                typename std::vector<int>::const_reference, // member type
                                                std::iter_reference_t<const_iterator_of_int_vector>, // iterator
                                                std::iter_reference_t<const_foreign_iterator>, // iterator2
                                                std::iter_reference_t<const_iterator_of_view>>; // range, no member

    using comp_list = seqan3::type_list<int const &,
                                        int const &,
                                        int const &,
                                        int const &,
                                        int>; // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, difference_type_)
{
    //TODO(h-2): add something that actually has a different difference_type
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using foreign_iterator = seqan3::detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = seqan3::type_list<std::ranges::range_difference_t<std::vector<int>>, // short
                                                typename std::vector<int>::difference_type, // member type
                                                std::ranges::range_difference_t<std::vector<int> const>, // const container
                                                std::iter_difference_t<iterator_of_int_vector>, // iterator
                                                std::iter_difference_t<foreign_iterator>, // iterator2
                                                std::ranges::range_difference_t<decltype(v)>>; // range, no member

    // views::ints' difference_type is not std::ptrdiff_t, but depends on the size.
    // For infinite views, like in our case, it's std::int_fast64_t (or std::int_fast32_t on 32bit).
    // On most platforms this is the same as std::size_t (`long int` on Linux and FreeBSD;
    // `long long int` on macOS and MSVC) BUT on some platforms like OpenBSD std::ptrdiff_t is
    // `long int`, but std::int64_t and std::int_fast64_t are `long long int`.
    // This detects word-size and selects the corresponding int_fast*_t which is always correct
    using view_int_diff_t = std::conditional_t<sizeof(std::ptrdiff_t) == 8,
                                               std::int_fast64_t,
                                               std::int_fast32_t>;

    using comp_list = seqan3::type_list<std::ptrdiff_t,
                                        std::ptrdiff_t,
                                        std::ptrdiff_t,
                                        std::ptrdiff_t,
                                        std::ptrdiff_t,
                                        view_int_diff_t>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, size_type_)
{
    //TODO(h-2): add something that actually has a different size_type
// iota is not a sized range, but take_exactly is
    auto v = std::views::iota(0) | seqan3::views::take_exactly(2);
    using type_list_example = seqan3::type_list<std::ranges::range_size_t<std::vector<int>>, // short
                                                typename std::vector<int>::size_type, // member type
                                                std::ranges::range_size_t<std::vector<int> const>, // const container
                                                std::ranges::range_size_t<decltype(v)>>; // range, no member

    using comp_list = seqan3::type_list<size_t,
                                        size_t,
                                        size_t,
                                        size_t>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, innermost_value_type_)
{
    using vector_of_int_vector = std::vector<std::vector<int>>;
    using iterator_of_int_vector = std::ranges::iterator_t<std::vector<int>>;
    using inner_value_type_of_const_iterator = seqan3::innermost_value_type_t<iterator_of_int_vector const>;
    using type_list_example = seqan3::type_list<typename seqan3::innermost_value_type<std::vector<int>>::type, // long
                                                seqan3::innermost_value_type_t<std::vector<int>>, // short
                                                seqan3::innermost_value_type_t<vector_of_int_vector>, // two-level
                                                seqan3::innermost_value_type_t<iterator_of_int_vector>, // iterator
                                                inner_value_type_of_const_iterator>; // const_iterator

    using comp_list = seqan3::type_list<int,
                                        int,
                                        int,
                                        int,
                                        int>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, dimension)
{
    EXPECT_EQ(1u, seqan3::dimension_v<std::vector<int>>);
    EXPECT_EQ(1u, seqan3::dimension_v<std::ranges::iterator_t<std::vector<int>>>);
    EXPECT_EQ(2u, seqan3::dimension_v<std::vector<std::vector<int>>>);
    EXPECT_EQ(2u, seqan3::dimension_v<std::ranges::iterator_t<std::vector<std::vector<int>>>>);
}

TEST(range_and_iterator, compatible)
{
    EXPECT_TRUE((seqan3::compatible<std::vector<int>,
                                    std::list<int>>));
    EXPECT_TRUE((seqan3::compatible<std::vector<int>,
                                    std::ranges::iterator_t<std::vector<int>>>));
    EXPECT_TRUE((seqan3::compatible<std::vector<int>,
                                    std::ranges::iterator_t<std::vector<int> const>>));
    EXPECT_TRUE((seqan3::compatible<std::list<std::vector<char>>,
                                    std::ranges::iterator_t<std::vector<std::string>>>));

    EXPECT_FALSE((seqan3::compatible<std::list<std::vector<char>>,
                                     std::string>));
    EXPECT_FALSE((seqan3::compatible<std::list<std::vector<char>>,
                                     std::ranges::iterator_t<std::string>>));
    EXPECT_FALSE((seqan3::compatible<std::list<int>,
                                     int>));
    EXPECT_FALSE((seqan3::compatible<std::vector<int>,
                                     std::string>));
}
