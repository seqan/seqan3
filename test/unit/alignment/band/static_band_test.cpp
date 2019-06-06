// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/band/static_band.hpp>

using namespace seqan3;

TEST(static_band, construction)
{
    static_band bs{lower_bound{-2}, upper_bound{2}};
    EXPECT_EQ(bs.lower_bound, -2);
    EXPECT_EQ(bs.upper_bound, 2);
}

TEST(static_band, wrong_boundary_args)
{
    EXPECT_THROW((static_band{lower_bound{3}, upper_bound{2}}), std::invalid_argument);
}
