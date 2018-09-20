// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <list>
#include <vector>

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/take_exactly.hpp>

#include <seqan3/core/metafunction/all.hpp>
<<<<<<< HEAD
=======
#include <seqan3/core/detail/reflection.hpp>
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
#include <seqan3/range/detail/random_access_iterator.hpp>

using namespace seqan3;

TEST(range_and_iterator, iterator_)
{
    EXPECT_TRUE((std::is_same_v<iterator_t<std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<iterator_t<std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));

    auto v = ranges::view::ints(1);
    EXPECT_TRUE((std::is_same_v<iterator_t<decltype(v)>,
                                decltype(ranges::begin(v))>));
    EXPECT_FALSE((std::is_same_v<iterator_t<decltype(v)>,
                                 decltype(ranges::end(v))>));
}

TEST(range_and_iterator, sentinel_)
{
    EXPECT_TRUE((std::is_same_v<sentinel_t<std::vector<int>>,
                                typename std::vector<int>::iterator>));
    EXPECT_TRUE((std::is_same_v<sentinel_t<std::vector<int> const>,
                                typename std::vector<int>::const_iterator>));
    EXPECT_TRUE((std::is_same_v<sentinel_t<std::vector<int>>,
                                iterator_t<std::vector<int>>>));

    auto v = ranges::view::ints(1);
    EXPECT_FALSE((std::is_same_v<sentinel_t<decltype(v)>,
                                decltype(ranges::begin(v))>));
    EXPECT_TRUE((std::is_same_v<sentinel_t<decltype(v)>,
                                decltype(ranges::end(v))>));
}

template <typename list1, typename list2, size_t pos = 0>
void expect_same_types()
{
<<<<<<< HEAD
    EXPECT_TRUE((std::is_same_v<meta::at_c<list1, pos>, meta::at_c<list2, pos>>));
    if constexpr (pos < list1::size() - 1)
        expect_same_types<list1, list2, pos + 1>();
    //NOTE(h-2): diagnostics could maybe be improved via get_display_name
=======
    constexpr bool val = std::is_same_v<meta::at_c<list1, pos>, meta::at_c<list2, pos>>;
    if constexpr (!val)
    {
        std::cerr << "\'" << detail::get_display_name_v<meta::at_c<list1, pos>>.string() << "\' not equal to \'"
                  << detail::get_display_name_v<meta::at_c<list2, pos>>.string() << "\' \n";
    }
    EXPECT_TRUE(val);

    if constexpr (pos < list1::size() - 1)
        expect_same_types<list1, list2, pos + 1>();
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
}

TEST(range_and_iterator, value_type_)
{
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = ranges::view::ints(1);
    using type_list = meta::list<value_type_t<std::vector<int>>,                    // short
                                 typename value_type<std::vector<int>>::type,       // long
                                 typename std::vector<int>::value_type,             // member type
                                 value_type_t<std::vector<int> const>,              // const container
                                 value_type_t<iterator_t<std::vector<int>>>,        // iterator
                                 value_type_t<foreign_iterator>,                    // iterator2
                                 value_type_t<decltype(v)>>;                        // range, no member
    using comp_list = meta::list<int,
                                 int,
                                 int,
                                 int,
                                 int,
                                 int,
                                 int>;

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, reference_)
{
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = ranges::view::ints(1);
    using type_list = meta::list<reference_t<std::vector<int>>,                    // short
                                 typename reference<std::vector<int>>::type,       // long
                                 typename std::vector<int>::reference,             // member type
                                 reference_t<std::vector<int> const>,              // const container
                                 reference_t<iterator_t<std::vector<int>>>,        // iterator
                                 reference_t<foreign_iterator>,                    // iterator2
                                 reference_t<decltype(v)>>;                        // range, no member

    using comp_list = meta::list<int &,
                                 int &,
                                 int &,
                                 int const &,   // container is const
                                 int &,
                                 int &,
                                 int>;          // view creates values

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, rvalue_reference_)
{
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = ranges::view::ints(1);
    using type_list = meta::list<rvalue_reference_t<std::vector<int>>,                    // short
                                 typename rvalue_reference<std::vector<int>>::type,       // long
// No types have member,yet:
//                                  typename std::vector<int>::rvalue_reference,             // member type
                                 rvalue_reference_t<std::vector<int> const>,              // const container
                                 rvalue_reference_t<iterator_t<std::vector<int>>>,        // iterator
                                 rvalue_reference_t<foreign_iterator>,                    // iterator2
                                 rvalue_reference_t<decltype(v)>>;                        // range, no member

    using comp_list = meta::list<int &&,
                                 int &&,
//                                  int &&,
                                 int const &&,  // container is const
                                 int &&,
                                 int &&,
                                 int>;          // view creates values

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, const_reference_)
{
//     using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = ranges::view::ints(1);
    using type_list = meta::list<const_reference_t<std::vector<int>>,                    // short
                                 typename const_reference<std::vector<int>>::type,       // long
                                 typename std::vector<int>::const_reference,             // member type
                                 const_reference_t<std::vector<int> const>,              // const container
// not defined on iterators
//                                  const_reference_t<iterator_t<std::vector<int>>>,        // iterator
//                                  const_reference_t<foreign_iterator>,                    // iterator2
                                 const_reference_t<decltype(v)>>;                        // range, no member

    using comp_list = meta::list<int const &,
                                 int const &,
                                 int const &,
                                 int const &,  // container is const
//                                  int const &,
//                                  int const &,
                                 int>;          // view creates values

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, difference_type_)
{
    //TODO(h-2): add something that actually has a different difference_type
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
    auto v = ranges::view::ints(1);
    using type_list = meta::list<difference_type_t<std::vector<int>>,                    // short
                                 typename difference_type<std::vector<int>>::type,       // long
                                 typename std::vector<int>::difference_type,             // member type
                                 difference_type_t<std::vector<int> const>,              // const container
                                 difference_type_t<iterator_t<std::vector<int>>>,        // iterator
                                 difference_type_t<foreign_iterator>,                    // iterator2
                                 difference_type_t<decltype(v)>>;                        // range, no member

<<<<<<< HEAD
    using comp_list = meta::list<std::make_signed_t<size_t>,
                                 std::make_signed_t<size_t>,
                                 std::make_signed_t<size_t>,
                                 std::make_signed_t<size_t>,
                                 std::make_signed_t<size_t>,
                                 std::make_signed_t<size_t>,
                                 std::make_signed_t<size_t>>;
=======
    // view::ints' difference_type is not std::ptrdiff_t, but depends on the size.
    // For infinite views, like in our case, it's std::int_fast64_t (or std::int_fast32_t on 32bit).
    // On most platforms this is the same as std::size_t (`long int` on Linux and FreeBSD;
    // `long long int` on macOS and MSVC) BUT on some platforms like OpenBSD std::ptrdiff_t is
    // `long int`, but std::int64_t and std::int_fast64_t are `long long int`.
    // This detects word-size and selects the corresponding int_fast*_t which is always correct
    using view_int_diff_t = std::conditional_t<sizeof(std::ptrdiff_t) == 8,
                                               std::int_fast64_t,
                                               std::int_fast32_t>;

    using comp_list = meta::list<std::ptrdiff_t,
                                 std::ptrdiff_t,
                                 std::ptrdiff_t,
                                 std::ptrdiff_t,
                                 std::ptrdiff_t,
                                 std::ptrdiff_t,
                                 view_int_diff_t>;
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, size_type_)
{
    //TODO(h-2): add something that actually has a different size_type
    using foreign_iterator = detail::random_access_iterator<std::vector<int>>;
// iota is not a sized range, but take_exactly is
    auto v = ranges::view::ints | ranges::view::take_exactly(2);
    using type_list = meta::list<size_type_t<std::vector<int>>,                    // short
                                 typename size_type<std::vector<int>>::type,       // long
                                 typename std::vector<int>::size_type,             // member type
                                 size_type_t<std::vector<int> const>,              // const container
                                 size_type_t<iterator_t<std::vector<int>>>,        // iterator
                                 size_type_t<foreign_iterator>,                    // iterator2
                                 size_type_t<decltype(v)>>;                        // range, no member

<<<<<<< HEAD
=======
    // see above
    using view_int_size_t = std::conditional_t<sizeof(std::size_t) == 8,
                                               std::uint_fast64_t,
                                               std::uint_fast32_t>;
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
    using comp_list = meta::list<size_t,
                                 size_t,
                                 size_t,
                                 size_t,
                                 size_t,
                                 size_t,
<<<<<<< HEAD
                                 size_t>;
=======
                                 view_int_size_t>;
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, innermost_value_type_)
{
    using type_list = meta::list<typename innermost_value_type<std::vector<int>>::type,         // long
                                 innermost_value_type_t<std::vector<int>>,                      // short
                                 innermost_value_type_t<std::vector<std::vector<int>>>,         // two-level
                                 innermost_value_type_t<iterator_t<std::vector<int>>>,          // iterator
                                 innermost_value_type_t<iterator_t<std::vector<int> const>>>;   // const_iterator

    using comp_list = meta::list<int,
                                 int,
                                 int,
                                 int,
                                 int>;

    expect_same_types<type_list, comp_list>();
}

TEST(range_and_iterator, dimension)
{
    EXPECT_EQ(dimension_v<std::vector<int>>,                            1u);
    EXPECT_EQ(dimension_v<iterator_t<std::vector<int>>>,                1u);
    EXPECT_EQ(dimension_v<std::vector<std::vector<int>>>,               2u);
    EXPECT_EQ(dimension_v<iterator_t<std::vector<std::vector<int>>>>,   2u);
}

TEST(range_and_iterator, compatible)
{
    EXPECT_TRUE((compatible_concept<std::vector<int>,
                                    std::list<int>>));
    EXPECT_TRUE((compatible_concept<std::vector<int>,
                                    iterator_t<std::vector<int>>>));
    EXPECT_TRUE((compatible_concept<std::vector<int>,
                                    iterator_t<std::vector<int> const>>));
    EXPECT_TRUE((compatible_concept<std::list<std::vector<char>>,
                                    iterator_t<std::vector<std::string>>>));

    EXPECT_FALSE((compatible_concept<std::list<std::vector<char>>,
                                     std::string>));
    EXPECT_FALSE((compatible_concept<std::list<std::vector<char>>,
                                     iterator_t<std::string>>));
    EXPECT_FALSE((compatible_concept<std::list<int>,
                                     int>));
    EXPECT_FALSE((compatible_concept<std::vector<int>,
                                     std::string>));
}
