// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <list>
#include <ranges>
#include <type_traits>
#include <vector>

#include <seqan3/utility/type_pack/traits.hpp>

// ----------------------------------------------------------------------------
// pack tests
// ----------------------------------------------------------------------------

TEST(pack_traits, size)
{
    EXPECT_EQ((seqan3::pack_traits::size<int, bool &, double const>), 3u);
}

TEST(pack_traits, count)
{
    EXPECT_EQ((seqan3::pack_traits::count<int>), 0u);
    EXPECT_EQ((seqan3::pack_traits::count<int, bool &, double const>), 0u);
    EXPECT_EQ((seqan3::pack_traits::count<int, bool &, int, double const, int>), 2u);
}

TEST(pack_traits, find)
{
    EXPECT_EQ((seqan3::pack_traits::find<int>), -1ll);
    EXPECT_EQ((seqan3::pack_traits::find<int, bool &, double const>), -1ll);
    EXPECT_EQ((seqan3::pack_traits::find<int, bool &, int, double const, int>), 1u);
}

TEST(pack_traits, find_if)
{
    EXPECT_EQ((seqan3::pack_traits::find_if<std::is_integral>), -1ll);
    EXPECT_EQ((seqan3::pack_traits::find_if<std::is_integral, float, double const>), -1ll);
    EXPECT_EQ((seqan3::pack_traits::find_if<std::is_integral, float, int, double const, long>), 1ll);
}

TEST(pack_traits, contains)
{
    EXPECT_EQ((seqan3::pack_traits::contains<int>), false);
    EXPECT_EQ((seqan3::pack_traits::contains<int, bool &, double const>), false);
    EXPECT_EQ((seqan3::pack_traits::contains<int, bool &, int, double const, int>), true);
}

TEST(pack_traits, at)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::at<2, int, bool &, double const, long, float>, double const>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::at<-2, int, bool &, double const, long, float>, long>));
}

TEST(pack_traits, front)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::front<int, bool &, double const, long, float>, int>));
}

TEST(pack_traits, back)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::back<int, bool &, double const, long, float>, float>));
}

TEST(pack_traits, drop_front)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::drop_front<int, bool &, double const, long, float>,
                                seqan3::type_list<bool &, double const, long, float>>));
}

TEST(pack_traits, transform)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::transform<std::ranges::range_value_t>, seqan3::type_list<>>));
    EXPECT_TRUE(
        (std::is_same_v<seqan3::pack_traits::transform<std::ranges::range_value_t, std::vector<int>, std::list<bool>>,
                        seqan3::type_list<int, bool>>));
    EXPECT_TRUE((std::is_same_v<
                 seqan3::pack_traits::transform<std::ranges::range_reference_t, std::vector<int>, std::list<bool>>,
                 seqan3::type_list<int &, bool &>>));
}

TEST(pack_traits, take)
{
    EXPECT_TRUE(
        (std::is_same_v<seqan3::pack_traits::take<0, int, bool &, double const, long, float>, seqan3::type_list<>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::take<3, int, bool &, double const, long, float>,
                                seqan3::type_list<int, bool &, double const>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::take<5, int, bool &, double const, long, float>,
                                seqan3::type_list<int, bool &, double const, long, float>>));
}

TEST(pack_traits, drop)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::drop<0, int, bool &, double const, long, float>,
                                seqan3::type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::drop<3, int, bool &, double const, long, float>,
                                seqan3::type_list<long, float>>));
    EXPECT_TRUE(
        (std::is_same_v<seqan3::pack_traits::drop<5, int, bool &, double const, long, float>, seqan3::type_list<>>));
}

TEST(pack_traits, take_last)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::take_last<0, int, bool &, double const, long, float>,
                                seqan3::type_list<>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::take_last<3, int, bool &, double const, long, float>,
                                seqan3::type_list<double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::take_last<5, int, bool &, double const, long, float>,
                                seqan3::type_list<int, bool &, double const, long, float>>));
}

TEST(pack_traits, drop_last)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::drop_last<0, int, bool &, double const, long, float>,
                                seqan3::type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::drop_last<3, int, bool &, double const, long, float>,
                                seqan3::type_list<int, bool &>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::drop_last<5, int, bool &, double const, long, float>,
                                seqan3::type_list<>>));
}

TEST(pack_traits, split_after)
{
    using split0 = seqan3::pack_traits::split_after<0, int, bool &, double const, long, float>;
    EXPECT_TRUE((std::is_same_v<typename split0::first_type, seqan3::type_list<>>));
    EXPECT_TRUE(
        (std::is_same_v<typename split0::second_type, seqan3::type_list<int, bool &, double const, long, float>>));

    using split3 = seqan3::pack_traits::split_after<3, int, bool &, double const, long, float>;
    EXPECT_TRUE((std::is_same_v<typename split3::first_type, seqan3::type_list<int, bool &, double const>>));
    EXPECT_TRUE((std::is_same_v<typename split3::second_type, seqan3::type_list<long, float>>));

    using split5 = seqan3::pack_traits::split_after<5, int, bool &, double const, long, float>;
    EXPECT_TRUE(
        (std::is_same_v<typename split5::first_type, seqan3::type_list<int, bool &, double const, long, float>>));
    EXPECT_TRUE((std::is_same_v<typename split5::second_type, seqan3::type_list<>>));
}

TEST(pack_traits, replace_at)
{
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::replace_at<double, 0, int, float, bool>,
                                seqan3::type_list<double, float, bool>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::replace_at<double, 1, int, float, bool>,
                                seqan3::type_list<int, double, bool>>));
    EXPECT_TRUE((std::is_same_v<seqan3::pack_traits::replace_at<double, 2, int, float, bool>,
                                seqan3::type_list<int, float, double>>));
}
