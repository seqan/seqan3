// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/ranges>
#include <string>

#include <range/v3/view/take.hpp>

TEST(ranges_test, combine_std_with_range_v3)
{
    std::string str{"foo"};
    auto take_first = str | std::views::take(5) | ranges::view::take(1);

    EXPECT_EQ(*std::ranges::begin(take_first), 'f');
}
