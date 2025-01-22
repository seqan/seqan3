// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/utility/detail/integer_traits.hpp>

TEST(int_types_test, min_viable_uint_t)
{
    using bool_1_t = seqan3::detail::min_viable_uint_t<0ull>;
    using bool_2_t = seqan3::detail::min_viable_uint_t<1ull>;
    using uint8_1_t = seqan3::detail::min_viable_uint_t<2ull>;
    using uint8_2_t = seqan3::detail::min_viable_uint_t<0xFFull>;
    using uint16_1_t = seqan3::detail::min_viable_uint_t<0x100ull>;
    using uint16_2_t = seqan3::detail::min_viable_uint_t<0xFF'FFull>;
    using uint32_1_t = seqan3::detail::min_viable_uint_t<0x1'00'00ull>;
    using uint32_2_t = seqan3::detail::min_viable_uint_t<0xFF'FF'FF'FFull>;
    using uint64_1_t = seqan3::detail::min_viable_uint_t<0x1'00'00'00'00ull>;
    using uint64_2_t = seqan3::detail::min_viable_uint_t<0xFF'FF'FF'FF'FF'FF'FF'FFull>;

    EXPECT_TRUE((std::is_same_v<bool_1_t, bool>));
    EXPECT_TRUE((std::is_same_v<bool_2_t, bool>));
    EXPECT_TRUE((std::is_same_v<uint8_1_t, uint8_t>));
    EXPECT_TRUE((std::is_same_v<uint8_2_t, uint8_t>));
    EXPECT_TRUE((std::is_same_v<uint16_1_t, uint16_t>));
    EXPECT_TRUE((std::is_same_v<uint16_2_t, uint16_t>));
    EXPECT_TRUE((std::is_same_v<uint32_1_t, uint32_t>));
    EXPECT_TRUE((std::is_same_v<uint32_2_t, uint32_t>));
    EXPECT_TRUE((std::is_same_v<uint64_1_t, uint64_t>));
    EXPECT_TRUE((std::is_same_v<uint64_2_t, uint64_t>));
}

TEST(int_types_test, min_viable_uint_v)
{
    auto bool_1_v = seqan3::detail::min_viable_uint_v<0ull>;
    auto bool_2_v = seqan3::detail::min_viable_uint_v<1ull>;
    auto uint8_1_v = seqan3::detail::min_viable_uint_v<2ull>;
    auto uint8_2_v = seqan3::detail::min_viable_uint_v<0xFFull>;
    auto uint16_1_v = seqan3::detail::min_viable_uint_v<0x100ull>;
    auto uint16_2_v = seqan3::detail::min_viable_uint_v<0xFF'FFull>;
    auto uint32_1_v = seqan3::detail::min_viable_uint_v<0x1'00'00ull>;
    auto uint32_2_v = seqan3::detail::min_viable_uint_v<0xFF'FF'FF'FFull>;
    auto uint64_1_v = seqan3::detail::min_viable_uint_v<0x1'00'00'00'00ull>;
    auto uint64_2_v = seqan3::detail::min_viable_uint_v<0xFF'FF'FF'FF'FF'FF'FF'FFull>;

    EXPECT_EQ(static_cast<uint64_t>(bool_1_v), 0ull);
    EXPECT_EQ(static_cast<uint64_t>(bool_2_v), 1ull);
    EXPECT_EQ(static_cast<uint64_t>(uint8_1_v), 2ull);
    EXPECT_EQ(static_cast<uint64_t>(uint8_2_v), 0xFFull);
    EXPECT_EQ(static_cast<uint64_t>(uint16_1_v), 0x100ull);
    EXPECT_EQ(static_cast<uint64_t>(uint16_2_v), 0xFF'FFull);
    EXPECT_EQ(static_cast<uint64_t>(uint32_1_v), 0x1'00'00ull);
    EXPECT_EQ(static_cast<uint64_t>(uint32_2_v), 0xFF'FF'FF'FFull);
    EXPECT_EQ(static_cast<uint64_t>(uint64_1_v), 0x1'00'00'00'00ull);
    EXPECT_EQ(static_cast<uint64_t>(uint64_2_v), 0xFF'FF'FF'FF'FF'FF'FF'FFull);

    EXPECT_TRUE((std::is_same_v<decltype(bool_1_v), bool>));
    EXPECT_TRUE((std::is_same_v<decltype(bool_2_v), bool>));
    EXPECT_TRUE((std::is_same_v<decltype(uint8_1_v), uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint8_2_v), uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint16_1_v), uint16_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint16_2_v), uint16_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint32_1_v), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint32_2_v), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint64_1_v), uint64_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint64_2_v), uint64_t>));
}
