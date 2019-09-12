// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/range/views/repeat.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/zip.hpp>

using namespace seqan3;

class zip_test : public ::testing::Test
{
protected:
    std::vector<int> vi{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<std::string> vs{"this", "is", "a", "test"};
    std::vector<std::string> const vc{"this", "is", "a", "test"};
    static constexpr auto vr = views::repeat('L');

    using tuple1_t = std::tuple<int, std::string, char>;
    using tuple2_t = std::tuple<int, std::string>;
    using result1_t = std::vector<tuple1_t>;
    using result2_t = std::vector<tuple2_t>;

    using rng_not_common_t = decltype(views::zip(vi, vs, vr));
    using rng_common_t = decltype(views::zip(vi, vs));
    using rng_const_t = decltype(views::zip(vi, vc));
    using view_const_t = decltype(views::zip(vi, vc)) const;
    using rng_view_const_t = decltype(views::zip(vi, vc)) const;

    result1_t expected1{{0, "this", 'L'}, {1, "is", 'L'}, {2, "a", 'L'}, {3, "test", 'L'}};
    result2_t expected2{{0, "this"}, {1, "is"}, {2, "a"}, {3, "test"}};

    rng_not_common_t v1;
    rng_common_t v2;

    virtual void SetUp()
    {
        v1 = views::zip(vi, vs, vr);
        v2 = views::zip(vi, vs);
    }
};

TEST_F(zip_test, concepts)
{
    EXPECT_TRUE((std::ranges::range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::range<rng_common_t>));
    EXPECT_TRUE((std::ranges::range<rng_const_t>));
    EXPECT_TRUE((std::ranges::range<view_const_t>));
    EXPECT_TRUE((std::ranges::range<rng_view_const_t>));

    EXPECT_TRUE((std::ranges::input_range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::input_range<rng_common_t>));
    EXPECT_TRUE((std::ranges::input_range<rng_const_t>));
    EXPECT_TRUE((std::ranges::input_range<view_const_t>));
    EXPECT_TRUE((std::ranges::input_range<rng_view_const_t>));

    EXPECT_TRUE((std::ranges::forward_range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::forward_range<rng_common_t>));
    EXPECT_TRUE((std::ranges::forward_range<rng_const_t>));
    EXPECT_TRUE((std::ranges::forward_range<view_const_t>));
    EXPECT_TRUE((std::ranges::forward_range<rng_view_const_t>));

    EXPECT_TRUE((std::ranges::bidirectional_range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::bidirectional_range<rng_common_t>));
    EXPECT_TRUE((std::ranges::bidirectional_range<rng_const_t>));
    EXPECT_TRUE((std::ranges::bidirectional_range<view_const_t>));
    EXPECT_TRUE((std::ranges::bidirectional_range<rng_view_const_t>));

    EXPECT_TRUE((std::ranges::random_access_range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::random_access_range<rng_common_t>));
    EXPECT_TRUE((std::ranges::random_access_range<rng_const_t>));
    EXPECT_TRUE((std::ranges::random_access_range<view_const_t>));
    EXPECT_TRUE((std::ranges::random_access_range<rng_view_const_t>));

    EXPECT_FALSE((std::ranges::contiguous_range<rng_not_common_t>));
    EXPECT_FALSE((std::ranges::contiguous_range<rng_common_t>));
    EXPECT_FALSE((std::ranges::contiguous_range<rng_const_t>));
    EXPECT_FALSE((std::ranges::contiguous_range<view_const_t>));
    EXPECT_FALSE((std::ranges::contiguous_range<rng_view_const_t>));

    EXPECT_TRUE((std::ranges::view<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::view<rng_common_t>));
    EXPECT_TRUE((std::ranges::view<rng_const_t>));
    EXPECT_FALSE((std::ranges::view<view_const_t>));
    EXPECT_FALSE((std::ranges::view<rng_view_const_t>));

    EXPECT_FALSE((std::ranges::sized_range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::sized_range<rng_common_t>));
    EXPECT_TRUE((std::ranges::sized_range<rng_const_t>));
    EXPECT_TRUE((std::ranges::sized_range<view_const_t>));
    EXPECT_TRUE((std::ranges::sized_range<rng_view_const_t>));

    EXPECT_FALSE((std::ranges::common_range<rng_not_common_t>));
    EXPECT_TRUE((std::ranges::common_range<rng_common_t>));
    EXPECT_TRUE((std::ranges::common_range<rng_const_t>));
    EXPECT_TRUE((std::ranges::common_range<view_const_t>));
    EXPECT_TRUE((std::ranges::common_range<rng_view_const_t>));

    EXPECT_TRUE((std::ranges::output_range<rng_not_common_t, std::tuple<int &, std::string &, char &>>));
    EXPECT_TRUE((std::ranges::output_range<rng_common_t, std::tuple<int &, std::string &>>));
    EXPECT_TRUE((std::ranges::output_range<rng_const_t, std::tuple<int &, std::string &>>));
    EXPECT_TRUE((std::ranges::output_range<view_const_t, std::tuple<int &, std::string &>>));
    EXPECT_TRUE((std::ranges::output_range<rng_view_const_t, std::tuple<int &, std::string &>>));
}

TEST_F(zip_test, access)
{
    EXPECT_TRUE(std::ranges::equal(v1, expected1));
    EXPECT_TRUE(std::ranges::equal(v2, expected2));
}

TEST_F(zip_test, combine)
{
    EXPECT_TRUE(std::ranges::equal(v1 | std::views::reverse, expected1 | std::views::reverse));
    EXPECT_TRUE(std::ranges::equal(v2 | std::views::reverse, expected2 | std::views::reverse));

    EXPECT_TRUE(std::ranges::equal(v1 | views::take(2), expected1 | views::take(2)));
    EXPECT_TRUE(std::ranges::equal(v1 | std::views::take(2), expected1 | std::views::take(2)));

    EXPECT_TRUE(std::ranges::equal(v2 | views::take(2), expected2 | views::take(2)));
    EXPECT_TRUE(std::ranges::equal(v2 | std::views::take(2), expected2 | std::views::take(2)));
}

TEST_F(zip_test, assign)
{
    *v1.begin() = std::tuple{9, "moo", 'P'};
    EXPECT_TRUE(std::ranges::equal(v1, result1_t{{9, "moo", 'P'}, {1, "is", 'P'}, {2, "a", 'P'}, {3, "test", 'P'}}));
}
