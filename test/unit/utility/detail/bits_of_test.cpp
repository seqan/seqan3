// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/utility/detail/bits_of.hpp>

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
