// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/std/algorithm>

TEST(general, move)
{
    std::vector v{0, 1, 2, 3, 4, 5, 6};
    std::vector<int> t;
    t.resize(v.size());

    std::ranges::move(v, t.begin());

    EXPECT_EQ(t, (std::vector{0, 1, 2, 3, 4, 5, 6}));
}

TEST(general, move_backward)
{
    std::vector v{0, 1, 2, 3, 4, 5, 6};
    std::vector<int> t;
    t.resize(v.size());

    std::ranges::move_backward(v, t.end());

    EXPECT_EQ(t, (std::vector{0, 1, 2, 3, 4, 5, 6}));
}
