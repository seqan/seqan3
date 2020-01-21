// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/std/algorithm>

TEST(general, find)
{
    std::vector v{0, 1, 2, 3, 4, 5, 6};

    auto res = std::ranges::find(v, 3);

    EXPECT_EQ(*res, 3);
}

TEST(general, find_if)
{
    std::vector v{0, 1, 2, 3, 4, 5, 6};

    auto res = std::ranges::find_if(v, [] (auto const c) { return c == 3; });

    EXPECT_EQ(*res, 3);
}

TEST(general, find_if_not)
{
    std::vector v{0, 1, 2, 3, 4, 5, 6};

    auto res = std::ranges::find_if_not(v, [] (auto const c) { return c != 3; });

    EXPECT_EQ(*res, 3);
}
