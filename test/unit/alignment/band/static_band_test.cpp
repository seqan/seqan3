// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/band/static_band.hpp>

TEST(static_band, construction)
{
    seqan3::static_band bs{seqan3::lower_bound{-2}, seqan3::upper_bound{2}};
    EXPECT_EQ(bs.lower_bound, -2);
    EXPECT_EQ(bs.upper_bound, 2);
}

TEST(static_band, wrong_boundary_args)
{
    EXPECT_THROW((seqan3::static_band{seqan3::lower_bound{3}, seqan3::upper_bound{2}}), std::invalid_argument);
}
