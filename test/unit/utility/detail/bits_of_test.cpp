// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/utility/detail/bits_of.hpp>

static constexpr size_t max_iterations = 1 << 15;

TEST(bits_of, bits_of)
{
    EXPECT_EQ(seqan3::detail::bits_of<int8_t>, 8);
    EXPECT_EQ(seqan3::detail::bits_of<uint8_t>, 8);
    EXPECT_EQ(seqan3::detail::bits_of<int16_t>, 16);
    EXPECT_EQ(seqan3::detail::bits_of<uint16_t>, 16);
    EXPECT_EQ(seqan3::detail::bits_of<int32_t>, 32);
    EXPECT_EQ(seqan3::detail::bits_of<uint32_t>, 32);
    EXPECT_EQ(seqan3::detail::bits_of<int64_t>, 64);
    EXPECT_EQ(seqan3::detail::bits_of<uint64_t>, 64);
}
