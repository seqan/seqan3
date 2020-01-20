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
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(range_and_iterator, iterator_)
{
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));

    auto v = std::views::iota(1);
    EXPECT_TRUE((std::is_same_v<std::ranges::iterator_t<decltype(v)>,
                                decltype(begin(v))>));
    EXPECT_FALSE((std::is_same_v<std::ranges::iterator_t<decltype(v)>,
                                 decltype(end(v))>));
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
                                decltype(begin(v))>));
    EXPECT_TRUE((std::is_same_v<std::ranges::sentinel_t<decltype(v)>,
                                decltype(end(v))>));
}

template <typename list1, typename list2, size_t pos = 0>
void expect_same_types()
{
    constexpr bool val = std::is_same_v<list_traits::at<pos, list1>, list_traits::at<pos, list2>>;
    if constexpr (!val)
    {
        std::cerr << "pos: " << pos << " \'" << detail::type_name_as_string<list_traits::at<pos, list1>>
                  << "\' not equal to \'" << detail::type_name_as_string<list_traits::at<pos, list2>> << "\' \n";
    }
    EXPECT_TRUE(val);

    if constexpr (pos < list1::size() - 1)
        expect_same_types<list1, list2, pos + 1>();
}

TEST(range_and_iterator, value_type_)
{
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = type_list<value_type_t<std::vector<int>>,                          // short
                                        typename value_type<std::vector<int>>::type,             // long
                                        typename std::vector<int>::value_type,                   // member type
                                        value_type_t<std::vector<int> const>,                    // const container
                                        value_type_t<std::ranges::iterator_t<std::vector<int>>>, // iterator
                                        value_type_t<foreign_iterator>,                          // iterator2
                                        value_type_t<decltype(v)>>;                              // range, no member
    using comp_list = type_list<int,
                                int,
                                int,
                                int,
                                int,
                                int,
                                int>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, reference_)
{
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = type_list<reference_t<std::vector<int>>,                          // short
                                        typename reference<std::vector<int>>::type,             // long
                                        typename std::vector<int>::reference,                   // member type
                                        reference_t<std::vector<int> const>,                    // const container
                                        reference_t<std::ranges::iterator_t<std::vector<int>>>, // iterator
                                        reference_t<foreign_iterator>,                          // iterator2
                                        reference_t<decltype(v)>>;                              // range, no member

    using comp_list = type_list<int &,
                                int &,
                                int &,
                                int const &,   // container is const
                                int &,
                                int &,
                                int>;          // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, rvalue_reference_)
{
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = type_list<rvalue_reference_t<std::vector<int>>,                          // short
                                        typename rvalue_reference<std::vector<int>>::type,             // long
// No types have member,yet:
//                                      typename std::vector<int>::rvalue_reference,                   // member type
                                        rvalue_reference_t<std::vector<int> const>,                    // const container
                                        rvalue_reference_t<std::ranges::iterator_t<std::vector<int>>>, // iterator
                                        rvalue_reference_t<foreign_iterator>,                          // iterator2
                                        rvalue_reference_t<decltype(v)>>;                              // range, no member

    using comp_list = type_list<int &&,
                                int &&,
//                              int &&,
                                int const &&,  // container is const
                                int &&,
                                int &&,
                                int>;          // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, const_reference_)
{
//     using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = type_list<const_reference_t<std::vector<int>>,                          // short
                                        typename const_reference<std::vector<int>>::type,             // long
                                        typename std::vector<int>::const_reference,                   // member type
                                        const_reference_t<std::vector<int> const>,                    // const container
// not defined on iterators
//                                      const_reference_t<std::ranges::iterator_t<std::vector<int>>>, // iterator
//                                      const_reference_t<foreign_iterator>,                          // iterator2
                                        const_reference_t<decltype(v)>>;                              // range, no member

    using comp_list = type_list<int const &,
                                int const &,
                                int const &,
                                int const &,  // container is const
//                              int const &,
//                              int const &,
                                int>;          // view creates values

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, difference_type_)
{
    //TODO(h-2): add something that actually has a different difference_type
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = std::views::iota(1);
    using type_list_example = type_list<difference_type_t<std::vector<int>>,                          // short
                                        typename difference_type<std::vector<int>>::type,             // long
                                        typename std::vector<int>::difference_type,                   // member type
                                        difference_type_t<std::vector<int> const>,                    // const container
                                        difference_type_t<std::ranges::iterator_t<std::vector<int>>>, // iterator
                                        difference_type_t<foreign_iterator>,                          // iterator2
                                        difference_type_t<decltype(v)>>;                              // range, no member

    // views::ints' difference_type is not std::ptrdiff_t, but depends on the size.
    // For infinite views, like in our case, it's std::int_fast64_t (or std::int_fast32_t on 32bit).
    // On most platforms this is the same as std::size_t (`long int` on Linux and FreeBSD;
    // `long long int` on macOS and MSVC) BUT on some platforms like OpenBSD std::ptrdiff_t is
    // `long int`, but std::int64_t and std::int_fast64_t are `long long int`.
    // This detects word-size and selects the corresponding int_fast*_t which is always correct
    using view_int_diff_t = std::conditional_t<sizeof(std::ptrdiff_t) == 8,
                                               std::int_fast64_t,
                                               std::int_fast32_t>;

    using comp_list = type_list<std::ptrdiff_t,
                                std::ptrdiff_t,
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
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
// iota is not a sized range, but take_exactly is
    auto v = std::views::iota(0) | views::take_exactly(2);
    using type_list_example = type_list<size_type_t<std::vector<int>>,                          // short
                                        typename size_type<std::vector<int>>::type,             // long
                                        typename std::vector<int>::size_type,                   // member type
                                        size_type_t<std::vector<int> const>,                    // const container
                                        size_type_t<std::ranges::iterator_t<std::vector<int>>>, // iterator
                                        size_type_t<foreign_iterator>,                          // iterator2
                                        size_type_t<decltype(v)>>;                              // range, no member

    using comp_list = type_list<size_t,
                                size_t,
                                size_t,
                                size_t,
                                size_t,
                                size_t,
                                size_t>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, innermost_value_type_)
{
    using type_list_example = type_list<typename innermost_value_type<std::vector<int>>::type,                    // long
                                        innermost_value_type_t<std::vector<int>>,                                 // short
                                        innermost_value_type_t<std::vector<std::vector<int>>>,                    // two-level
                                        innermost_value_type_t<std::ranges::iterator_t<std::vector<int>>>,        // iterator
                                        innermost_value_type_t<std::ranges::iterator_t<std::vector<int> const>>>; // const_iterator

    using comp_list = type_list<int,
                                int,
                                int,
                                int,
                                int>;

    expect_same_types<type_list_example, comp_list>();
}

TEST(range_and_iterator, dimension)
{
    EXPECT_EQ(dimension_v<std::vector<int>>,                            1u);
    EXPECT_EQ(dimension_v<std::ranges::iterator_t<std::vector<int>>>,                1u);
    EXPECT_EQ(dimension_v<std::vector<std::vector<int>>>,               2u);
    EXPECT_EQ(dimension_v<std::ranges::iterator_t<std::vector<std::vector<int>>>>,   2u);
}

TEST(range_and_iterator, compatible)
{
    EXPECT_TRUE((compatible<std::vector<int>,
                                    std::list<int>>));
    EXPECT_TRUE((compatible<std::vector<int>,
                                    std::ranges::iterator_t<std::vector<int>>>));
    EXPECT_TRUE((compatible<std::vector<int>,
                                    std::ranges::iterator_t<std::vector<int> const>>));
    EXPECT_TRUE((compatible<std::list<std::vector<char>>,
                                    std::ranges::iterator_t<std::vector<std::string>>>));

    EXPECT_FALSE((compatible<std::list<std::vector<char>>,
                                     std::string>));
    EXPECT_FALSE((compatible<std::list<std::vector<char>>,
                                     std::ranges::iterator_t<std::string>>));
    EXPECT_FALSE((compatible<std::list<int>,
                                     int>));
    EXPECT_FALSE((compatible<std::vector<int>,
                                     std::string>));
}
