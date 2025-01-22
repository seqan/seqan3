// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/utility/detail/to_little_endian.hpp>

TEST(to_little_endian, byte)
{
    uint8_t val = 0x01;

    EXPECT_EQ(seqan3::detail::to_little_endian(val), 0x01);
}

TEST(to_little_endian, word)
{
    uint16_t val = 0x01'02; // 258
    uint16_t res = seqan3::detail::to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x02'); // LSB
    EXPECT_EQ(*(res_p + 1), '\x01'); // MSB
}

TEST(to_little_endian, double_word)
{
    uint32_t val = 0x01'02'03'04; // 16.909.060
    uint32_t res = seqan3::detail::to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x04');
    EXPECT_EQ(*(res_p + 1), '\x03');
    EXPECT_EQ(*(res_p + 2), '\x02');
    EXPECT_EQ(*(res_p + 3), '\x01');
}

TEST(to_little_endian, quad_word)
{
    uint64_t val = 0x01'02'03'04'05'06'07'08;
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
