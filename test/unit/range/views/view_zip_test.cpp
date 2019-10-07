// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/range/views/repeat.hpp>
#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/zip.hpp>

#include <range/v3/view/drop_exactly.hpp>



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
    // using t [[maybe_unused]] = decltype(views::zip(std::declval<std::vector<int> &>(), std::declval<std::vector<int> &>()));
    // using test [[maybe_unused]] = decltype(views::single_pass_input(views::zip(std::declval<std::vector<int> &>(), views::repeat('L'))));
    // using resource_t [[maybe_unused]] = std::span<std::tuple<int, int>>;
    // using calamitas [[maybe_unused]] = decltype(views::single_pass_input(views::zip(std::declval<resource_t>(), std::views::iota(1))));
    // using minimal_calamitas [[maybe_unused]] = decltype(views::single_pass_input(views::zip(views::repeat('l'))));
    // using not_calamitas [[maybe_unused]] = decltype(views::single_pass_input(views::zip(std::declval<resource_t>(), std::declval<std::vector<int> &>())));

    auto v9 = views::single_pass_input(views::zip(std::views::iota(0)));
    using T = decltype(v9);
    using it_t = std::ranges::iterator_t<T>;
    auto it = v9.begin();
    auto it2 = v9.end();
    (void) it;
    (void) it2;
    static_assert(std::same_as<decltype(*std::declval<it_t>()), std::ranges::iter_reference_t<it_t>>);
    static_assert(std::same_as<decltype(ranges::_::iter_move(std::declval<it_t>())), std::ranges::iter_rvalue_reference_t<it_t>>);
    // using minimal_calamitas [[maybe_unused]] = decltype(views::single_pass_input(views::zip(std::views::iota(0))));
    // requires_<same_as<decltype (ranges::_::iter_move(i)), ranges::iter_rvalue_reference_t<I> > >)

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
