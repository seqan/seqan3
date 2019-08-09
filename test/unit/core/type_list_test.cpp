// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/core/type_list/traits.hpp>

using namespace seqan3;

TEST(type_list, basic)
{
    using t = type_list<int, char, double>;

    EXPECT_TRUE((std::is_same_v<meta::at_c<t, 1>, char>));
}

// ----------------------------------------------------------------------------
// pack tests
// ----------------------------------------------------------------------------

TEST(pack_traits, size)
{
    EXPECT_EQ((pack_traits::size<int, bool &, double const>),
               3u);
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

// ----------------------------------------------------------------------------
// list tests
// ----------------------------------------------------------------------------

TEST(list_traits, size)
{
    EXPECT_EQ((list_traits::size<type_list<int, bool &, double const>>),
               3u);
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
