// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/bit_manipulation.hpp>

static constexpr size_t max_iterations = 1 << 15;

TEST(bit_manipulation, sizeof_bits)
{
    EXPECT_EQ(seqan3::detail::sizeof_bits<int8_t>, 8);
    EXPECT_EQ(seqan3::detail::sizeof_bits<uint8_t>, 8);
    EXPECT_EQ(seqan3::detail::sizeof_bits<int16_t>, 16);
    EXPECT_EQ(seqan3::detail::sizeof_bits<uint16_t>, 16);
    EXPECT_EQ(seqan3::detail::sizeof_bits<int32_t>, 32);
    EXPECT_EQ(seqan3::detail::sizeof_bits<uint32_t>, 32);
    EXPECT_EQ(seqan3::detail::sizeof_bits<int64_t>, 64);
    EXPECT_EQ(seqan3::detail::sizeof_bits<uint64_t>, 64);
}

TEST(bit_manipulation, is_power_of_two)
{
    constexpr bool is_power_of_two0 = seqan3::detail::is_power_of_two(0);
    constexpr bool is_power_of_two1 = seqan3::detail::is_power_of_two(1);
    constexpr bool is_power_of_two2 = seqan3::detail::is_power_of_two(2);
    constexpr bool is_power_of_two3 = seqan3::detail::is_power_of_two(3);
    EXPECT_FALSE(is_power_of_two0);
    EXPECT_TRUE(is_power_of_two1);
    EXPECT_TRUE(is_power_of_two2);
    EXPECT_FALSE(is_power_of_two3);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_TRUE(seqan3::detail::is_power_of_two(power_of_two));

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_FALSE(seqan3::detail::is_power_of_two(i)) << i << " should not be a power of two.";
        }
    }
}

TEST(bit_manipulation, next_power_of_two)
{
    constexpr size_t next_power_of_two0 = seqan3::detail::next_power_of_two(0);
    constexpr size_t next_power_of_two1 = seqan3::detail::next_power_of_two(1);
    constexpr size_t next_power_of_two2 = seqan3::detail::next_power_of_two(2);
    constexpr size_t next_power_of_two3 = seqan3::detail::next_power_of_two(3);
    EXPECT_EQ(next_power_of_two0, 1u);
    EXPECT_EQ(next_power_of_two1, 1u);
    EXPECT_EQ(next_power_of_two2, 2u);
    EXPECT_EQ(next_power_of_two3, 4u);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_EQ(seqan3::detail::next_power_of_two(power_of_two), power_of_two);

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_EQ(seqan3::detail::next_power_of_two(i), next_power) << "The next power of two of " << i
                                                                        << " should be " << next_power;
        }
    }
}

using unsigned_types = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;

template <typename type>
class unsigned_operations : public ::testing::Test
{};

TYPED_TEST_SUITE(unsigned_operations, unsigned_types, );

TYPED_TEST(unsigned_operations, most_significant_bit_set)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero  = seqan3::detail::most_significant_bit_set<unsigned_t>(0b0001);
    constexpr size_t one1  = seqan3::detail::most_significant_bit_set<unsigned_t>(0b0010);
    constexpr size_t one2  = seqan3::detail::most_significant_bit_set<unsigned_t>(0b0011);
    constexpr size_t two1  = seqan3::detail::most_significant_bit_set<unsigned_t>(0b0101);
    constexpr size_t two2  = seqan3::detail::most_significant_bit_set<unsigned_t>(0b0111);
    constexpr size_t seven = seqan3::detail::most_significant_bit_set<unsigned_t>(0b10010010);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(one1, 1u);
    EXPECT_EQ(one2, 1u);
    EXPECT_EQ(two1, 2u);
    EXPECT_EQ(two2, 2u);
    EXPECT_EQ(seven, 7u);

    for (uint8_t position = 0; position < seqan3::detail::sizeof_bits<unsigned_t>; ++position)
    {
        unsigned_t start = unsigned_t{1u} << position;
        unsigned_t end = start << 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(sdsl::bits::hi(n), position) << "[SDSL] The position of the msb of " << n << " should be "
                                                   << position;
            EXPECT_EQ(seqan3::detail::most_significant_bit_set(n), position) << "The position of the msb of " << n
                                                                             << " should be " << position;
        }
    }
}

TYPED_TEST(unsigned_operations, count_leading_zeros)
{
    using unsigned_t = TypeParam;
    constexpr size_t t1 = seqan3::detail::count_leading_zeros<unsigned_t>(0b0001);
    constexpr size_t t2 = seqan3::detail::count_leading_zeros<unsigned_t>(0b0101);
    constexpr size_t t3 = seqan3::detail::count_leading_zeros<unsigned_t>(0b0010);
    constexpr size_t t4 = seqan3::detail::count_leading_zeros<unsigned_t>(0b0110);
    constexpr size_t t5 = seqan3::detail::count_leading_zeros<unsigned_t>(0b0100);
    constexpr size_t t6 = seqan3::detail::count_leading_zeros<unsigned_t>(0b10100000);
    EXPECT_EQ(t1, seqan3::detail::sizeof_bits<unsigned_t> - 1u);
    EXPECT_EQ(t2, seqan3::detail::sizeof_bits<unsigned_t> - 3u);
    EXPECT_EQ(t3, seqan3::detail::sizeof_bits<unsigned_t> - 2u);
    EXPECT_EQ(t4, seqan3::detail::sizeof_bits<unsigned_t> - 3u);
    EXPECT_EQ(t5, seqan3::detail::sizeof_bits<unsigned_t> - 3u);
    EXPECT_EQ(t6, seqan3::detail::sizeof_bits<unsigned_t> - 8u);

    for (uint8_t cnt = 0; cnt < seqan3::detail::sizeof_bits<unsigned_t>; ++cnt)
    {
        unsigned_t start = std::numeric_limits<unsigned_t>::max() >> cnt;
        unsigned_t end = start >> 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(seqan3::detail::sizeof_bits<unsigned_t> - sdsl::bits::hi(n) - 1, n) << "[SDSL] n " << n
                                                                                          << " should have " << cnt
                                                                                          << " leading zeros.";
            EXPECT_EQ(seqan3::detail::count_leading_zeros(n), cnt) << "n " << n << " should have " << cnt
                                                                   << " leading zeros.";
        }
    }
}

TYPED_TEST(unsigned_operations, count_trailing_zeros)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero  = seqan3::detail::count_trailing_zeros<unsigned_t>(0b0001);
    constexpr size_t zero2 = seqan3::detail::count_trailing_zeros<unsigned_t>(0b0101);
    constexpr size_t one1  = seqan3::detail::count_trailing_zeros<unsigned_t>(0b0010);
    constexpr size_t one2  = seqan3::detail::count_trailing_zeros<unsigned_t>(0b0110);
    constexpr size_t two   = seqan3::detail::count_trailing_zeros<unsigned_t>(0b0100);
    constexpr size_t five  = seqan3::detail::count_trailing_zeros<unsigned_t>(0b10100000);
    EXPECT_EQ(zero,  0u);
    EXPECT_EQ(zero2, 0u);
    EXPECT_EQ(one1,  1u);
    EXPECT_EQ(one2,  1u);
    EXPECT_EQ(two,   2u);
    EXPECT_EQ(five,  5u);

    for (uint8_t cnt = 0; cnt < seqan3::detail::sizeof_bits<unsigned_t>; ++cnt)
    {
        unsigned_t start = std::numeric_limits<unsigned_t>::max() << cnt;
        unsigned_t end = start << 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(sdsl::bits::lo(n), cnt) << "[SDSL] n " << n << " should have " << cnt << " trailing zeros.";
            EXPECT_EQ(seqan3::detail::count_trailing_zeros(n), cnt) << "n " << n << " should have " << cnt
                                                                    << " trailing zeros.";
        }
    }
}

// https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
template <std::unsigned_integral unsigned_t>
unsigned_t permute_bits(unsigned_t v)
{
    if (v & (unsigned_t{1u} << (seqan3::detail::sizeof_bits<unsigned_t> - 1)))
        return v;

    unsigned_t t = v | (v - 1);
    unsigned_t w = (t + 1) | (((~t & -~t) - 1) >> (seqan3::detail::count_trailing_zeros(v) + 1));
    return w;
}

TYPED_TEST(unsigned_operations, popcount)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero  = seqan3::detail::popcount<unsigned_t>(0b0000);
    constexpr size_t one   = seqan3::detail::popcount<unsigned_t>(0b0100);
    constexpr size_t two   = seqan3::detail::popcount<unsigned_t>(0b1100);
    constexpr size_t three = seqan3::detail::popcount<unsigned_t>(0b1110);
    constexpr size_t four  = seqan3::detail::popcount<unsigned_t>(0b1111);
    constexpr size_t five  = seqan3::detail::popcount<unsigned_t>(0b10011011);
    EXPECT_EQ(zero,  0u);
    EXPECT_EQ(one,   1u);
    EXPECT_EQ(two,   2u);
    EXPECT_EQ(three, 3u);
    EXPECT_EQ(four,  4u);
    EXPECT_EQ(five,  5u);

    for (uint8_t position = 0; position < seqan3::detail::sizeof_bits<unsigned_t>; ++position)
    {
        unsigned_t start = std::numeric_limits<unsigned_t>::max() >> position;
        auto sizeof_bits_of_unsigned_t = seqan3::detail::sizeof_bits<unsigned_t>;

        EXPECT_EQ(seqan3::detail::popcount(start),
                  sizeof_bits_of_unsigned_t - position) << "The pocount of " << start << " should be "
                                                        << sizeof_bits_of_unsigned_t - position;
        for (unsigned_t n = permute_bits(start), k = 0u;
             n > start && k < max_iterations;
             start = n, n = permute_bits(start), ++k)
        {
            EXPECT_EQ(static_cast<uint8_t>(sdsl::bits::cnt(n)),
                      sizeof_bits_of_unsigned_t - position) << "[SDSL] The pocount of " << n << " should be "
                                                            << sizeof_bits_of_unsigned_t - position;
            EXPECT_EQ(seqan3::detail::popcount(n),
                      sizeof_bits_of_unsigned_t - position) << "The pocount of " << n << " should be "
                                                            << sizeof_bits_of_unsigned_t - position;
        }
    }
}

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
