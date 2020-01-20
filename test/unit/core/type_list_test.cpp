// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <list>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_list/all.hpp>

using namespace seqan3;

TEST(type_list, basic)
{
    using t = type_list<int, char, double>;

    EXPECT_TRUE((std::is_same_v<list_traits::at<1, t>, char>));
}

// ----------------------------------------------------------------------------
// pack tests
// ----------------------------------------------------------------------------

TEST(pack_traits, size)
{
    EXPECT_EQ((pack_traits::size<int, bool &, double const>),
               3u);
}

TEST(pack_traits, count)
{
    EXPECT_EQ((pack_traits::count<int>),
               0u);
    EXPECT_EQ((pack_traits::count<int, bool &, double const>),
               0u);
    EXPECT_EQ((pack_traits::count<int, bool &, int, double const, int>),
               2u);
}

TEST(pack_traits, find)
{
    EXPECT_EQ((pack_traits::find<int>),
               -1ll);
    EXPECT_EQ((pack_traits::find<int, bool &, double const>),
               -1ll);
    EXPECT_EQ((pack_traits::find<int, bool &, int, double const, int>),
               1u);
}

TEST(pack_traits, find_if)
{
    EXPECT_EQ((pack_traits::find_if<std::is_integral>),
               -1ll);
    EXPECT_EQ((pack_traits::find_if<std::is_integral, float, double const>),
               -1ll);
    EXPECT_EQ((pack_traits::find_if<std::is_integral, float, int, double const, long>),
               1ll);
}

TEST(pack_traits, contains)
{
    EXPECT_EQ((pack_traits::contains<int>),
               false);
    EXPECT_EQ((pack_traits::contains<int, bool &, double const>),
               false);
    EXPECT_EQ((pack_traits::contains<int, bool &, int, double const, int>),
               true);
}

TEST(pack_traits, at)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::at<2, int, bool &, double const, long, float>,
                                double const>));
    EXPECT_TRUE((std::is_same_v<pack_traits::at<-2, int, bool &, double const, long, float>,
                                long>));
}

TEST(pack_traits, front)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::front<int, bool &, double const, long, float>,
                                int>));
}

TEST(pack_traits, back)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::back<int, bool &, double const, long, float>,
                                float>));
}

TEST(pack_traits, drop_front)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::drop_front<int, bool &, double const, long, float>,
                                type_list<bool &, double const, long, float>>));
}

TEST(pack_traits, transform)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::transform<value_type_t>,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::transform<value_type_t, std::vector<int>, std::list<bool>>,
                                type_list<int, bool>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::transform<reference_t, std::vector<int>, std::list<bool>>,
                                type_list<int &, bool &>>));
}

TEST(pack_traits, take)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::take<0, int, bool &, double const, long, float>,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::take<3, int, bool &, double const, long, float>,
                                type_list<int, bool &, double const>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::take<5, int, bool &, double const, long, float>,
                                type_list<int, bool &, double const, long, float>>));
}

TEST(pack_traits, drop)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::drop<0, int, bool &, double const, long, float>,
                                type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::drop<3, int, bool &, double const, long, float>,
                                type_list<long, float>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::drop<5, int, bool &, double const, long, float>,
                                type_list<>>));
}

TEST(pack_traits, take_last)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::take_last<0, int, bool &, double const, long, float>,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::take_last<3, int, bool &, double const, long, float>,
                                type_list<double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::take_last<5, int, bool &, double const, long, float>,
                                type_list<int, bool &, double const, long, float>>));
}

TEST(pack_traits, drop_last)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::drop_last<0, int, bool &, double const, long, float>,
                                type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::drop_last<3, int, bool &, double const, long, float>,
                                type_list<int, bool &>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::drop_last<5, int, bool &, double const, long, float>,
                                type_list<>>));
}

TEST(pack_traits, split_after)
{
    using split0 = pack_traits::split_after<0, int, bool &, double const, long, float>;
    EXPECT_TRUE((std::is_same_v<typename split0::first_type,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<typename split0::second_type,
                                type_list<int, bool &, double const, long, float>>));

    using split3 = pack_traits::split_after<3, int, bool &, double const, long, float>;
    EXPECT_TRUE((std::is_same_v<typename split3::first_type,
                                type_list<int, bool &, double const>>));
    EXPECT_TRUE((std::is_same_v<typename split3::second_type,
                                type_list<long, float>>));


    using split5 = pack_traits::split_after<5, int, bool &, double const, long, float>;
    EXPECT_TRUE((std::is_same_v<typename split5::first_type,
                                type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<typename split5::second_type,
                                type_list<>>));
}

TEST(pack_traits, replace_at)
{
    EXPECT_TRUE((std::is_same_v<pack_traits::replace_at<double, 0, int, float, bool>,
                                type_list<double, float, bool>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::replace_at<double, 1, int, float, bool>,
                                type_list<int, double, bool>>));
    EXPECT_TRUE((std::is_same_v<pack_traits::replace_at<double, 2, int, float, bool>,
                                type_list<int, float, double>>));
}

// ----------------------------------------------------------------------------
// list tests
// ----------------------------------------------------------------------------

TEST(list_traits, size)
{
    EXPECT_EQ((list_traits::size<type_list<int, bool &, double const>>),
               3u);
}

TEST(list_traits, count)
{
    EXPECT_EQ((list_traits::count<int, type_list<>>),
               0u);
    EXPECT_EQ((list_traits::count<int, type_list<bool &, double const>>),
               0u);
    EXPECT_EQ((list_traits::count<int, type_list<bool &, int, double const, int>>),
               2u);
}

TEST(list_traits, find)
{
    EXPECT_EQ((list_traits::find<int, type_list<>>),
               -1ll);
    EXPECT_EQ((list_traits::find<int, type_list<bool &, double const>>),
               -1ll);
    EXPECT_EQ((list_traits::find<int, type_list<bool &, int, double const, int>>),
               1u);
}

TEST(list_traits, find_if)
{
    EXPECT_EQ((list_traits::find_if<std::is_integral, type_list<>>),
               -1ll);
    EXPECT_EQ((list_traits::find_if<std::is_integral, type_list<float, double const>>),
               -1ll);
    EXPECT_EQ((list_traits::find_if<std::is_integral, type_list<float, int, double const, long>>),
               1ll);
}

TEST(list_traits, contains)
{
    EXPECT_EQ((list_traits::contains<int, type_list<>>),
               false);
    EXPECT_EQ((list_traits::contains<int, type_list<bool &, double const>>),
               false);
    EXPECT_EQ((list_traits::contains<int, type_list<bool &, int, double const, int>>),
               true);
}

TEST(list_traits, at)
{
    EXPECT_TRUE((std::is_same_v<list_traits::at<2, type_list<int, bool &, double const, long, float>>,
                                double const>));
    EXPECT_TRUE((std::is_same_v<list_traits::at<-2, type_list<int, bool &, double const, long, float>>,
                                long>));
}

TEST(list_traits, front)
{
    EXPECT_TRUE((std::is_same_v<list_traits::front<type_list<int, bool &, double const, long, float>>,
                                int>));
}

TEST(list_traits, back)
{
    EXPECT_TRUE((std::is_same_v<list_traits::back<type_list<int, bool &, double const, long, float>>,
                                float>));
}

TEST(list_traits, concat)
{
    EXPECT_TRUE((std::is_same_v<list_traits::concat<type_list<int, bool &, double const>, type_list<long, float>>,
                                type_list<int, bool &, double const, long, float>>));

    EXPECT_TRUE((std::is_same_v<list_traits::concat<type_list<int, bool &, double const>,
                                                    type_list<long, float>,
                                                    type_list<>,
                                                    type_list<long &>>,
                                type_list<int, bool &, double const, long, float, long &>>));
}

TEST(list_traits, drop_front)
{
    EXPECT_TRUE((std::is_same_v<list_traits::drop_front<type_list<int, bool &, double const, long, float>>,
                                type_list<bool &, double const, long, float>>));
}

TEST(list_traits, take)
{
    EXPECT_TRUE((std::is_same_v<list_traits::take<0, type_list<int, bool &, double const, long, float>>,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<list_traits::take<3, type_list<int, bool &, double const, long, float>>,
                                type_list<int, bool &, double const>>));
    EXPECT_TRUE((std::is_same_v<list_traits::take<5, type_list<int, bool &, double const, long, float>>,
                                type_list<int, bool &, double const, long, float>>));
}

TEST(list_traits, drop)
{
    EXPECT_TRUE((std::is_same_v<list_traits::drop<0, type_list<int, bool &, double const, long, float>>,
                                type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<list_traits::drop<3, type_list<int, bool &, double const, long, float>>,
                                type_list<long, float>>));
    EXPECT_TRUE((std::is_same_v<list_traits::drop<5, type_list<int, bool &, double const, long, float>>,
                                type_list<>>));
}

TEST(list_traits, take_last)
{
    EXPECT_TRUE((std::is_same_v<list_traits::take_last<0, type_list<int, bool &, double const, long, float>>,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<list_traits::take_last<3, type_list<int, bool &, double const, long, float>>,
                                type_list<double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<list_traits::take_last<5, type_list<int, bool &, double const, long, float>>,
                                type_list<int, bool &, double const, long, float>>));
}

TEST(list_traits, drop_last)
{
    EXPECT_TRUE((std::is_same_v<list_traits::drop_last<0, type_list<int, bool &, double const, long, float>>,
                                type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<list_traits::drop_last<3, type_list<int, bool &, double const, long, float>>,
                                type_list<int, bool &>>));
    EXPECT_TRUE((std::is_same_v<list_traits::drop_last<5, type_list<int, bool &, double const, long, float>>,
                                type_list<>>));
}

TEST(list_traits, split_after)
{
    using split0 = list_traits::split_after<0, type_list<int, bool &, double const, long, float>>;
    EXPECT_TRUE((std::is_same_v<typename split0::first_type,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<typename split0::second_type,
                                type_list<int, bool &, double const, long, float>>));

    using split3 = list_traits::split_after<3, type_list<int, bool &, double const, long, float>>;
    EXPECT_TRUE((std::is_same_v<typename split3::first_type,
                                type_list<int, bool &, double const>>));
    EXPECT_TRUE((std::is_same_v<typename split3::second_type,
                                type_list<long, float>>));


    using split5 = list_traits::split_after<5, type_list<int, bool &, double const, long, float>>;
    EXPECT_TRUE((std::is_same_v<typename split5::first_type,
                                type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<typename split5::second_type,
                                type_list<>>));
}

TEST(list_traits, transform)
{
    EXPECT_TRUE((std::is_same_v<list_traits::transform<value_type_t, type_list<>>,
                                type_list<>>));
    EXPECT_TRUE((std::is_same_v<list_traits::transform<value_type_t, type_list<std::vector<int>, std::list<bool>>>,
                                type_list<int, bool>>));
    EXPECT_TRUE((std::is_same_v<list_traits::transform<reference_t, type_list<std::vector<int>, std::list<bool>>>,
                                type_list<int &, bool &>>));
}

TEST(list_traits, replace_at)
{
    EXPECT_TRUE((std::is_same_v<list_traits::replace_at<double, 0, type_list<int, float, bool>>,
                                type_list<double, float, bool>>));
    EXPECT_TRUE((std::is_same_v<list_traits::replace_at<double, 1, type_list<int, float, bool>>,
                                type_list<int, double, bool>>));
    EXPECT_TRUE((std::is_same_v<list_traits::replace_at<double, 2, type_list<int, float, bool>>,
                                type_list<int, float, double>>));
}
