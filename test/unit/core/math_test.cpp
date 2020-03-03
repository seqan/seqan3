// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/math.hpp>

TEST(pow, unsigned_base)
{
    EXPECT_EQ(0u, seqan3::pow(0u, 2u));
    EXPECT_EQ(1u, seqan3::pow(2u, 0u));
    EXPECT_EQ(8u, seqan3::pow(2u, 3u));
    EXPECT_EQ(std::numeric_limits<uint64_t>::max(), seqan3::pow(std::numeric_limits<uint64_t>::max(), 1u));
}

TEST(pow, signed_base)
{
    EXPECT_EQ(0, seqan3::pow(0, 2u));
    EXPECT_EQ(1, seqan3::pow(2, 0u));
    EXPECT_EQ(8, seqan3::pow(2, 3u));
    EXPECT_EQ(-8, seqan3::pow(-2, 3u));
    EXPECT_EQ(std::numeric_limits<int64_t>::max(), seqan3::pow(std::numeric_limits<int64_t>::max(), 1u));
    EXPECT_EQ(std::numeric_limits<int64_t>::min(), seqan3::pow(std::numeric_limits<int64_t>::min(), 1u));
}

TEST(pow, std)
{
    EXPECT_EQ(0.0, seqan3::pow(0u, 2));
    EXPECT_EQ(1.0, seqan3::pow(2, 0));
    EXPECT_EQ(27.0, seqan3::pow(3.0, 3u));
    EXPECT_EQ(-8.0, seqan3::pow(-2.0, 3));
}

#ifndef NDEBUG
TEST(pow, overflow)
{
    EXPECT_THROW(seqan3::pow(2u, 64u), std::overflow_error);
    EXPECT_THROW(seqan3::pow(2, 63u), std::overflow_error);
}

TEST(pow, underflow)
{
    EXPECT_THROW(seqan3::pow(-3, 50u), std::underflow_error);
}
#endif
