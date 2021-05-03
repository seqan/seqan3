// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <list>
#include <seqan3/std/ranges>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/type_list/traits.hpp>

TEST(list_traits, size)
{
    EXPECT_EQ((seqan3::list_traits::size<seqan3::type_list<int, bool &, double const>>),
               3u);
}

TEST(list_traits, count)
{
    EXPECT_EQ((seqan3::list_traits::count<int, seqan3::type_list<>>),
               0u);
    EXPECT_EQ((seqan3::list_traits::count<int, seqan3::type_list<bool &, double const>>),
               0u);
    EXPECT_EQ((seqan3::list_traits::count<int, seqan3::type_list<bool &, int, double const, int>>),
               2u);
}

TEST(list_traits, find)
{
    EXPECT_EQ((seqan3::list_traits::find<int, seqan3::type_list<>>),
               -1ll);
    EXPECT_EQ((seqan3::list_traits::find<int, seqan3::type_list<bool &, double const>>),
               -1ll);
    EXPECT_EQ((seqan3::list_traits::find<int, seqan3::type_list<bool &, int, double const, int>>),
               1u);
}

TEST(list_traits, find_if)
{
    EXPECT_EQ((seqan3::list_traits::find_if<std::is_integral, seqan3::type_list<>>),
               -1ll);
    EXPECT_EQ((seqan3::list_traits::find_if<std::is_integral, seqan3::type_list<float, double const>>),
               -1ll);
    EXPECT_EQ((seqan3::list_traits::find_if<std::is_integral, seqan3::type_list<float, int, double const, long>>),
               1ll);
}

TEST(list_traits, contains)
{
    EXPECT_EQ((seqan3::list_traits::contains<int, seqan3::type_list<>>),
               false);
    EXPECT_EQ((seqan3::list_traits::contains<int, seqan3::type_list<bool &, double const>>),
               false);
    EXPECT_EQ((seqan3::list_traits::contains<int, seqan3::type_list<bool &, int, double const, int>>),
               true);
}

TEST(list_traits, at)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE((seqan3::list_traits::at<2, test_types_list>), double const);
    EXPECT_SAME_TYPE((seqan3::list_traits::at<-2, test_types_list>), long);
}

TEST(list_traits, front)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE(seqan3::list_traits::front<test_types_list>, int);
}

TEST(list_traits, back)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE(seqan3::list_traits::back<test_types_list>, float);
}

TEST(list_traits, concat)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE((seqan3::list_traits::concat<seqan3::type_list<int, bool &, double const>,
                                                  seqan3::type_list<long, float>>),
                     test_types_list);

    EXPECT_SAME_TYPE((seqan3::list_traits::concat<seqan3::type_list<int, bool &, double const>,
                                                  seqan3::type_list<long, float>,
                                                  seqan3::type_list<>,
                                                  seqan3::type_list<long &>>),
                     (seqan3::type_list<int, bool &, double const, long, float, long &>));
}

TEST(list_traits, drop_front)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE(seqan3::list_traits::drop_front<test_types_list>,
                     (seqan3::type_list<bool &, double const, long, float>));
}

TEST(list_traits, take)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE((seqan3::list_traits::take<0, test_types_list>), seqan3::type_list<>);
    EXPECT_SAME_TYPE((seqan3::list_traits::take<3, test_types_list>), (seqan3::type_list<int, bool &, double const>));
    EXPECT_SAME_TYPE((seqan3::list_traits::take<5, test_types_list>), test_types_list);
}

TEST(list_traits, drop)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE((seqan3::list_traits::drop<0, test_types_list>), test_types_list);
    EXPECT_SAME_TYPE((seqan3::list_traits::drop<3, test_types_list>), (seqan3::type_list<long, float>));
    EXPECT_SAME_TYPE((seqan3::list_traits::drop<5, test_types_list>), seqan3::type_list<>);
}

TEST(list_traits, take_last)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE((seqan3::list_traits::take_last<0, test_types_list>), seqan3::type_list<>);
    EXPECT_SAME_TYPE((seqan3::list_traits::take_last<3, test_types_list>),
                     (seqan3::type_list<double const, long, float>));
    EXPECT_SAME_TYPE((seqan3::list_traits::take_last<5, test_types_list>), test_types_list);
}

TEST(list_traits, drop_last)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    EXPECT_SAME_TYPE((seqan3::list_traits::drop_last<0, test_types_list>), test_types_list);
    EXPECT_SAME_TYPE((seqan3::list_traits::drop_last<3, test_types_list>), (seqan3::type_list<int, bool &>));
    EXPECT_SAME_TYPE((seqan3::list_traits::drop_last<5, test_types_list>), seqan3::type_list<>);
}

TEST(list_traits, split_after)
{
    using test_types_list = seqan3::type_list<int, bool &, double const, long, float>;
    using split0 = seqan3::list_traits::split_after<0, test_types_list>;
    EXPECT_SAME_TYPE(typename split0::first_type, seqan3::type_list<>);
    EXPECT_SAME_TYPE(typename split0::second_type, test_types_list);

    using split3 = seqan3::list_traits::split_after<3, test_types_list>;
    EXPECT_SAME_TYPE(typename split3::first_type, (seqan3::type_list<int, bool &, double const>));
    EXPECT_SAME_TYPE(typename split3::second_type, (seqan3::type_list<long, float>));


    using split5 = seqan3::list_traits::split_after<5, test_types_list>;
    EXPECT_SAME_TYPE(typename split5::first_type, test_types_list);
    EXPECT_SAME_TYPE(typename split5::second_type, seqan3::type_list<>);
}

TEST(list_traits, transform)
{
    EXPECT_SAME_TYPE((seqan3::list_traits::transform<std::ranges::range_value_t,
                                                     seqan3::type_list<>>),
                     seqan3::type_list<>);
    EXPECT_SAME_TYPE((seqan3::list_traits::transform<std::ranges::range_value_t,
                                                     seqan3::type_list<std::vector<int>, std::list<bool>>>),
                     (seqan3::type_list<int, bool>));
    EXPECT_SAME_TYPE((seqan3::list_traits::transform<std::ranges::range_reference_t,
                                                     seqan3::type_list<std::vector<int>, std::list<bool>>>),
                     (seqan3::type_list<int &, bool &>));
}

TEST(list_traits, replace_at)
{
    EXPECT_SAME_TYPE((seqan3::list_traits::replace_at<double, 0, seqan3::type_list<int, float, bool>>),
                     (seqan3::type_list<double, float, bool>));
    EXPECT_SAME_TYPE((seqan3::list_traits::replace_at<double, 1, seqan3::type_list<int, float, bool>>),
                     (seqan3::type_list<int, double, bool>));
    EXPECT_SAME_TYPE((seqan3::list_traits::replace_at<double, 2, seqan3::type_list<int, float, bool>>),
                     (seqan3::type_list<int, float, double>));
}

TEST(list_traits_detail, reverse)
{
    auto reverse = [] (auto && type_list)
    {
        return seqan3::list_traits::detail::reverse(type_list);
    };

    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<>{})), (seqan3::type_list<>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<float>{})), (seqan3::type_list<float>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<float, double, char, short>{})),
                     (seqan3::type_list<short, char, double, float>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<int>{})), (seqan3::type_list<int>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<int, int, int, int>{})),
                     (seqan3::type_list<int, int, int, int>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<float, int>{})), (seqan3::type_list<int, float>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<int, float>{})), (seqan3::type_list<float, int>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<int, float, int, int, double, int, char, short, int>{})),
                     (seqan3::type_list<int, short, char, int, double, int, int, float, int>));
    EXPECT_SAME_TYPE(decltype(reverse(seqan3::type_list<int, int, int, int, float, float, float>{})),
                     (seqan3::type_list<float, float, float, int, int, int, int>));
}
