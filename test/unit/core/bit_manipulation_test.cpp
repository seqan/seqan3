// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/bit_manipulation.hpp>

using namespace seqan3::detail;

static constexpr size_t max_iterations = 1 << 15;

TEST(bit_manipulation, sizeof_bits)
{
    EXPECT_EQ(sizeof_bits<int8_t>, 8);
    EXPECT_EQ(sizeof_bits<uint8_t>, 8);
    EXPECT_EQ(sizeof_bits<int16_t>, 16);
    EXPECT_EQ(sizeof_bits<uint16_t>, 16);
    EXPECT_EQ(sizeof_bits<int32_t>, 32);
    EXPECT_EQ(sizeof_bits<uint32_t>, 32);
    EXPECT_EQ(sizeof_bits<int64_t>, 64);
    EXPECT_EQ(sizeof_bits<uint64_t>, 64);
}

TEST(bit_manipulation, is_power_of_two)
{
    constexpr bool is_power_of_two0 = is_power_of_two(0);
    constexpr bool is_power_of_two1 = is_power_of_two(1);
    constexpr bool is_power_of_two2 = is_power_of_two(2);
    constexpr bool is_power_of_two3 = is_power_of_two(3);
    EXPECT_FALSE(is_power_of_two0);
    EXPECT_TRUE(is_power_of_two1);
    EXPECT_TRUE(is_power_of_two2);
    EXPECT_FALSE(is_power_of_two3);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_TRUE(is_power_of_two(power_of_two));

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_FALSE(is_power_of_two(i)) << i << " should not be a power of two.";
        }
    }
}

TEST(bit_manipulation, next_power_of_two)
{
    constexpr size_t next_power_of_two0 = next_power_of_two(0);
    constexpr size_t next_power_of_two1 = next_power_of_two(1);
    constexpr size_t next_power_of_two2 = next_power_of_two(2);
    constexpr size_t next_power_of_two3 = next_power_of_two(3);
    EXPECT_EQ(next_power_of_two0, 1u);
    EXPECT_EQ(next_power_of_two1, 1u);
    EXPECT_EQ(next_power_of_two2, 2u);
    EXPECT_EQ(next_power_of_two3, 4u);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_EQ(next_power_of_two(power_of_two), power_of_two);

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_EQ(next_power_of_two(i), next_power) << "The next power of two of " << i << " should be " << next_power;
        }
    }
}

using unsigned_types = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;

template <typename type>
class unsigned_operations : public ::testing::Test
{};

TYPED_TEST_CASE(unsigned_operations, unsigned_types);

TYPED_TEST(unsigned_operations, bit_scan_reverse)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero = bit_scan_reverse<unsigned_t>(0b0001);
    constexpr size_t one1 = bit_scan_reverse<unsigned_t>(0b0010);
    constexpr size_t one2 = bit_scan_reverse<unsigned_t>(0b0011);
    constexpr size_t two1 = bit_scan_reverse<unsigned_t>(0b0101);
    constexpr size_t two2 = bit_scan_reverse<unsigned_t>(0b0111);
    constexpr size_t seven = bit_scan_reverse<unsigned_t>(0b10010010);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(one1, 1u);
    EXPECT_EQ(one2, 1u);
    EXPECT_EQ(two1, 2u);
    EXPECT_EQ(two2, 2u);
    EXPECT_EQ(seven, 7u);

    for (uint8_t position = 0; position < 8u * sizeof(unsigned_t); ++position)
    {
        unsigned_t start = unsigned_t{1u} << position;
        unsigned_t end = start << 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(bit_scan_reverse(n), position) << "The position of the msb of " << n << " should be " << position;
        }
    }
}

TEST(to_little_endian, byte)
{
    uint8_t val = 0x01;

    EXPECT_EQ(to_little_endian(val), 0x01);
}

TEST(to_little_endian, word)
{
    uint16_t val = 0x0102;  // 258
    uint16_t res = to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x02'); // LSB
    EXPECT_EQ(*(res_p + 1), '\x01'); // MSB
}

TEST(to_little_endian, double_word)
{
    uint32_t val = 0x01020304; // 16.909.060
    uint32_t res = to_little_endian(val);
    char * res_p = reinterpret_cast<char *>(&res);

    EXPECT_EQ(*(res_p + 0), '\x04');
    EXPECT_EQ(*(res_p + 1), '\x03');
    EXPECT_EQ(*(res_p + 2), '\x02');
    EXPECT_EQ(*(res_p + 3), '\x01');
}

TEST(to_little_endian, quad_word)
{
    uint64_t val = 0x0102030405060708;
    uint64_t res = to_little_endian(val);
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
