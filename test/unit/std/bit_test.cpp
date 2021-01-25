// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/bit>
#include <seqan3/std/concepts>

#include <gtest/gtest.h>

#include <sdsl/bits.hpp>

#include <seqan3/utility/detail/bits_of.hpp>

static constexpr size_t max_iterations = 1 << 15;

TEST(bit, has_single_bit)
{
    constexpr bool is_power_of_two0 = std::has_single_bit(0u);
    constexpr bool is_power_of_two1 = std::has_single_bit(1u);
    constexpr bool is_power_of_two2 = std::has_single_bit(2u);
    constexpr bool is_power_of_two3 = std::has_single_bit(3u);
    EXPECT_FALSE(is_power_of_two0);
    EXPECT_TRUE(is_power_of_two1);
    EXPECT_TRUE(is_power_of_two2);
    EXPECT_FALSE(is_power_of_two3);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_TRUE(std::has_single_bit(power_of_two));

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_FALSE(std::has_single_bit(i)) << i << " should not be a power of two.";
        }
    }
}

TEST(bit, bit_ceil)
{
    constexpr size_t next_power_of_two0 = std::bit_ceil(0u);
    constexpr size_t next_power_of_two1 = std::bit_ceil(1u);
    constexpr size_t next_power_of_two2 = std::bit_ceil(2u);
    constexpr size_t next_power_of_two3 = std::bit_ceil(3u);
    EXPECT_EQ(next_power_of_two0, 1u);
    EXPECT_EQ(next_power_of_two1, 1u);
    EXPECT_EQ(next_power_of_two2, 2u);
    EXPECT_EQ(next_power_of_two3, 4u);

    for (size_t power_of_two = 1; power_of_two <= (size_t{1u} << 31); power_of_two <<= 1)
    {
        EXPECT_EQ(std::bit_ceil(power_of_two), power_of_two);

        size_t next_power = (power_of_two << 1);
        for (size_t i = power_of_two + 1, k = 0; i < next_power && k < max_iterations; ++i, ++k)
        {
            EXPECT_EQ(std::bit_ceil(i), next_power) << "The next power of two of " << i << " should be " << next_power;
        }
    }
}

using unsigned_types = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;

template <typename type>
class unsigned_operations : public ::testing::Test
{};

TYPED_TEST_SUITE(unsigned_operations, unsigned_types, );

TYPED_TEST(unsigned_operations, bit_width)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero = std::bit_width<unsigned_t>(0b0000);
    constexpr size_t one = std::bit_width<unsigned_t>(0b0001);
    constexpr size_t two1 = std::bit_width<unsigned_t>(0b0010);
    constexpr size_t two2 = std::bit_width<unsigned_t>(0b0011);
    constexpr size_t three1 = std::bit_width<unsigned_t>(0b0101);
    constexpr size_t three2 = std::bit_width<unsigned_t>(0b0111);
    constexpr size_t eight = std::bit_width<unsigned_t>(0b10010010);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(one, 1u);
    EXPECT_EQ(two1, 2u);
    EXPECT_EQ(two2, 2u);
    EXPECT_EQ(three1, 3u);
    EXPECT_EQ(three2, 3u);
    EXPECT_EQ(eight, 8u);

    for (uint8_t position = 0; position < seqan3::detail::bits_of<unsigned_t>; ++position)
    {
        unsigned_t start = unsigned_t{1u} << position;
        unsigned_t end = start << 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(sdsl::bits::hi(n), position) << "[SDSL] The position of the msb of " << n << " should be "
                                                   << position;
            EXPECT_EQ(std::bit_width(n), position + 1u) << "The position of the msb of " << n << " should be "
                                                        << position;
        }
    }
}

TYPED_TEST(unsigned_operations, countl_zero)
{
    using unsigned_t = TypeParam;
    constexpr size_t t0 = std::countl_zero<unsigned_t>(0b0000u);
    constexpr size_t t1 = std::countl_zero<unsigned_t>(0b0001u);
    constexpr size_t t2 = std::countl_zero<unsigned_t>(0b0101u);
    constexpr size_t t3 = std::countl_zero<unsigned_t>(0b0010u);
    constexpr size_t t4 = std::countl_zero<unsigned_t>(0b0110u);
    constexpr size_t t5 = std::countl_zero<unsigned_t>(0b0100u);
    constexpr size_t t6 = std::countl_zero<unsigned_t>(0b10100000u);
    EXPECT_EQ(t0, seqan3::detail::bits_of<unsigned_t>);
    EXPECT_EQ(t1, seqan3::detail::bits_of<unsigned_t> - 1u);
    EXPECT_EQ(t2, seqan3::detail::bits_of<unsigned_t> - 3u);
    EXPECT_EQ(t3, seqan3::detail::bits_of<unsigned_t> - 2u);
    EXPECT_EQ(t4, seqan3::detail::bits_of<unsigned_t> - 3u);
    EXPECT_EQ(t5, seqan3::detail::bits_of<unsigned_t> - 3u);
    EXPECT_EQ(t6, seqan3::detail::bits_of<unsigned_t> - 8u);

    for (uint8_t cnt = 0; cnt < seqan3::detail::bits_of<unsigned_t>; ++cnt)
    {
        unsigned_t start = std::numeric_limits<unsigned_t>::max() >> cnt;
        unsigned_t end = start >> 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(seqan3::detail::bits_of<unsigned_t> - sdsl::bits::hi(n) - 1, n) << "[SDSL] n " << n
                                                                                      << " should have " << cnt
                                                                                      << " leading zeros.";
            EXPECT_EQ(std::countl_zero(n), cnt) << "n " << n << " should have " << cnt
                                                             << " leading zeros.";
            EXPECT_EQ(std::countl_zero(n), cnt);
        }
    }
}

TYPED_TEST(unsigned_operations, countr_zero)
{
    using unsigned_t = TypeParam;
    constexpr size_t bits_of = std::countr_zero<unsigned_t>(0b0000);
    constexpr size_t zero = std::countr_zero<unsigned_t>(0b0001);
    constexpr size_t zero2 = std::countr_zero<unsigned_t>(0b0101);
    constexpr size_t one1 = std::countr_zero<unsigned_t>(0b0010);
    constexpr size_t one2 = std::countr_zero<unsigned_t>(0b0110);
    constexpr size_t two = std::countr_zero<unsigned_t>(0b0100);
    constexpr size_t five = std::countr_zero<unsigned_t>(0b10100000);
    EXPECT_EQ(bits_of, seqan3::detail::bits_of<unsigned_t>);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(zero2, 0u);
    EXPECT_EQ(one1, 1u);
    EXPECT_EQ(one2, 1u);
    EXPECT_EQ(two, 2u);
    EXPECT_EQ(five, 5u);

    for (uint8_t cnt = 0; cnt < seqan3::detail::bits_of<unsigned_t>; ++cnt)
    {
        unsigned_t start = std::numeric_limits<unsigned_t>::max() << cnt;
        unsigned_t end = start << 1u;
        for (unsigned_t n = start, k = 0u; n < end && k < max_iterations; ++n, ++k)
        {
            EXPECT_EQ(sdsl::bits::lo(n), cnt) << "[SDSL] n " << n << " should have " << cnt << " trailing zeros.";
            EXPECT_EQ(std::countr_zero(n), cnt) << "n " << n << " should have " << cnt << " trailing zeros.";
        }
    }
}

// https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
template <std::unsigned_integral unsigned_t>
unsigned_t permute_bits(unsigned_t v)
{
    if (v & (unsigned_t{1u} << (seqan3::detail::bits_of<unsigned_t> - 1)))
        return v;

    unsigned_t t = v | (v - 1);
    unsigned_t w = (t + 1) | (((~t & -~t) - 1) >> (std::countr_zero(v) + 1));
    return w;
}

TYPED_TEST(unsigned_operations, popcount)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero = std::popcount<unsigned_t>(0b0000);
    constexpr size_t one = std::popcount<unsigned_t>(0b0100);
    constexpr size_t two = std::popcount<unsigned_t>(0b1100);
    constexpr size_t three = std::popcount<unsigned_t>(0b1110);
    constexpr size_t four = std::popcount<unsigned_t>(0b1111);
    constexpr size_t five = std::popcount<unsigned_t>(0b10011011);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(one, 1u);
    EXPECT_EQ(two, 2u);
    EXPECT_EQ(three, 3u);
    EXPECT_EQ(four, 4u);
    EXPECT_EQ(five, 5u);

    for (uint8_t position = 0; position < seqan3::detail::bits_of<unsigned_t>; ++position)
    {
        unsigned_t start = std::numeric_limits<unsigned_t>::max() >> position;
        auto sizeof_bits_of_unsigned_t = seqan3::detail::bits_of<unsigned_t>;

        EXPECT_EQ(std::popcount(start),
                  sizeof_bits_of_unsigned_t - position) << "The pocount of " << start << " should be "
                                                        << sizeof_bits_of_unsigned_t - position;
        for (unsigned_t n = permute_bits(start), k = 0u;
             n > start && k < max_iterations;
             start = n, n = permute_bits(start), ++k)
        {
            EXPECT_EQ(static_cast<uint8_t>(sdsl::bits::cnt(n)),
                      sizeof_bits_of_unsigned_t - position) << "[SDSL] The pocount of " << n << " should be "
                                                            << sizeof_bits_of_unsigned_t - position;
            EXPECT_EQ(std::popcount(n),
                      sizeof_bits_of_unsigned_t - position) << "The pocount of " << n << " should be "
                                                            << sizeof_bits_of_unsigned_t - position;
            EXPECT_EQ(std::popcount(n), sizeof_bits_of_unsigned_t - position);
        }
    }
}
