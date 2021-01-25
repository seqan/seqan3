// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/utility/detail/to_little_endian.hpp>

TEST(to_little_endian, byte)
{
    uint8_t val = 0x01;

    EXPECT_EQ(seqan3::detail::to_little_endian(val), 0x01);
}

TEST(to_little_endian, word)
{
    uint16_t val = 0x0102;  // 258
    uint16_t res = seqan3::detail::to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x02'); // LSB
    EXPECT_EQ(*(res_p + 1), '\x01'); // MSB
}

TEST(to_little_endian, double_word)
{
    uint32_t val = 0x01020304; // 16.909.060
    uint32_t res = seqan3::detail::to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x04');
    EXPECT_EQ(*(res_p + 1), '\x03');
    EXPECT_EQ(*(res_p + 2), '\x02');
    EXPECT_EQ(*(res_p + 3), '\x01');
}

TEST(to_little_endian, quad_word)
{
    uint64_t val = 0x0102030405060708;
    uint64_t res = seqan3::detail::to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x08');
    EXPECT_EQ(*(res_p + 1), '\x07');
    EXPECT_EQ(*(res_p + 2), '\x06');
    EXPECT_EQ(*(res_p + 3), '\x05');
    EXPECT_EQ(*(res_p + 4), '\x04');
    EXPECT_EQ(*(res_p + 5), '\x03');
    EXPECT_EQ(*(res_p + 6), '\x02');
    EXPECT_EQ(*(res_p + 7), '\x01');
}
