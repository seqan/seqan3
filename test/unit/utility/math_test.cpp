// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/utility/detail/bits_of.hpp>
#include <seqan3/utility/math.hpp>

static constexpr size_t max_iterations = 1 << 15;

using unsigned_types = ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t>;

template <typename type>
class unsigned_operations : public ::testing::Test
{};

TYPED_TEST_SUITE(unsigned_operations, unsigned_types, );

TYPED_TEST(unsigned_operations, floor_log2)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero = seqan3::detail::floor_log2<unsigned_t>(0b0001);
    constexpr size_t one1 = seqan3::detail::floor_log2<unsigned_t>(0b0010);
    constexpr size_t one2 = seqan3::detail::floor_log2<unsigned_t>(0b0011);
    constexpr size_t two1 = seqan3::detail::floor_log2<unsigned_t>(0b0101);
    constexpr size_t two2 = seqan3::detail::floor_log2<unsigned_t>(0b0111);
    constexpr size_t seven = seqan3::detail::floor_log2<unsigned_t>(0b1001'0010);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(one1, 1u);
    EXPECT_EQ(one2, 1u);
    EXPECT_EQ(two1, 2u);
    EXPECT_EQ(two2, 2u);
    EXPECT_EQ(seven, 7u);

    for (uint8_t log2_value = 0; log2_value < seqan3::detail::bits_of<unsigned_t>; ++log2_value)
    {
        unsigned_t start = unsigned_t{1u} << log2_value;
        unsigned_t end = static_cast<unsigned_t>(std::min<size_t>(start, max_iterations) + start);
        for (unsigned_t n = start; n < end; ++n)
        {
            EXPECT_EQ(seqan3::detail::floor_log2(n), log2_value);
            EXPECT_EQ(std::floor(std::log2(n)), log2_value) << "If this fails this might be a floating point rounding "
                                                            << "error on your machine";
        }
    }
}

TYPED_TEST(unsigned_operations, ceil_log2)
{
    using unsigned_t = TypeParam;
    constexpr size_t zero = seqan3::detail::ceil_log2<unsigned_t>(0b0001);
    constexpr size_t one = seqan3::detail::ceil_log2<unsigned_t>(0b0010);
    constexpr size_t two = seqan3::detail::ceil_log2<unsigned_t>(0b0011);
    constexpr size_t three1 = seqan3::detail::ceil_log2<unsigned_t>(0b0101);
    constexpr size_t three2 = seqan3::detail::ceil_log2<unsigned_t>(0b0111);
    constexpr size_t eight = seqan3::detail::ceil_log2<unsigned_t>(0b1001'0010);
    EXPECT_EQ(zero, 0u);
    EXPECT_EQ(one, 1u);
    EXPECT_EQ(two, 2u);
    EXPECT_EQ(three1, 3u);
    EXPECT_EQ(three2, 3u);
    EXPECT_EQ(eight, 8u);

    for (uint8_t log2_value = 0; log2_value < seqan3::detail::bits_of<unsigned_t>; ++log2_value)
    {
        unsigned_t start = unsigned_t{1u} << log2_value;
        unsigned_t end = static_cast<unsigned_t>(std::min<size_t>(start, max_iterations) + start);
        EXPECT_EQ(seqan3::detail::ceil_log2(start), log2_value);
        EXPECT_EQ(std::ceil(std::log2(start)), log2_value) << "ceil_log2 of " << start << " should be " << log2_value
                                                           << "; If this fails this might be a floating point rounding "
                                                           << "error on your machine.";

        for (unsigned_t n = start + 1u; n < end; ++n)
        {
            EXPECT_EQ(seqan3::detail::ceil_log2(n), log2_value + 1u);

            if constexpr (seqan3::detail::bits_of<unsigned_t> <= 32u) // known to fail for 64bit unsigned integers
            {
                EXPECT_EQ(std::ceil(std::log2(n)), log2_value + 1u)
                    << "ceil_log2 of " << start << " should be " << log2_value
                    << "; If this fails this might be a floating point"
                    << "rounding error on your machine.";
            }
        }
    }
}

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
