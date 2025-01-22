// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

// make sure that including the std header does not produce any errors
// see https://github.com/seqan/seqan3/issues/2352
#include <charconv>
#include <cmath>
#include <concepts>
#include <iostream>
#include <limits>
#include <seqan3/std/charconv>

// =============================================================================
// std::from_chars for integral types
// =============================================================================

template <typename T>
class integral_from_char_test : public ::testing::Test
{};

using integral_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t>;

TYPED_TEST_SUITE(integral_from_char_test, integral_types, );

TYPED_TEST(integral_from_char_test, postive_number)
{
    TypeParam value{42};

    {
        std::vector<char> const str{'1', '2', '3'};
        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(value, TypeParam{123});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
        EXPECT_EQ(res.ec, std::errc{});
    }

    {
        std::vector<char> const str{'0', '2', '3'};
        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(value, TypeParam{23});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
        EXPECT_EQ(res.ec, std::errc{});
    }

    // Read only up to a certain point
    {
        std::vector<char> const str{'0', '2', '3', '4', '5', '6'};
        auto res = std::from_chars(&str[0], &str[0] + 3, value);

        EXPECT_EQ(value, TypeParam{23});
        EXPECT_EQ(res.ptr, &str[0] + 3);
        EXPECT_EQ(res.ec, std::errc{});
    }
}

TYPED_TEST(integral_from_char_test, negative_number)
{
    TypeParam value{42};
    std::vector<char> const str{'-', '1', '2', '3'};
    auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

    if constexpr (std::unsigned_integral<TypeParam>)
    {
        EXPECT_EQ(res.ptr, &str[0]);
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
        EXPECT_EQ(value, TypeParam{42});
    }
    else
    {
        EXPECT_EQ(value, TypeParam{-123});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
        EXPECT_EQ(res.ec, std::errc{});
    }
}

TYPED_TEST(integral_from_char_test, overflow_error)
{
    TypeParam value{42};
    std::vector<char> const str{'1', '2', '3', '0', '0', '0', '0', '0', '0', '0', '0',
                                '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};

    auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

    EXPECT_EQ(res.ptr, &str[0] + str.size());
    EXPECT_EQ(res.ec, std::errc::result_out_of_range);
    EXPECT_EQ(value, TypeParam{42});
}

TYPED_TEST(integral_from_char_test, partial_parsing)
{
    TypeParam value{42};

    { // interleaved char
        std::vector<char> const str{'1', 'a', '3'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(res.ptr, &str[1]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{1});
    }

    value = 42; // reset
    {           // trailing char
        std::vector<char> const str{'1', '2', 'a'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(res.ptr, &str[2]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{12});
    }

    value = 42; // reset
    {           // float
        std::vector<char> const str{'1', '.', '3'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(res.ptr, &str[1]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{1});
    }

    { // hexadecimal 0x prefix is not recognized
        std::vector<char> const str{'0', 'x', '3', 'f'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value, 16);

        EXPECT_EQ(res.ptr, &str[1]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{0});
    }

    { // hexadecimal 0X prefix is not recognized
        std::vector<char> const str{'0', 'X', '3', 'f'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value, 16);

        EXPECT_EQ(res.ptr, &str[1]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{0});
    }
}

TYPED_TEST(integral_from_char_test, invalid_argument_error)
{
    TypeParam value{42};

    { // leading char
        std::vector<char> const str{'a', '1', '3'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(res.ptr, &str[0]);
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
        EXPECT_EQ(value, TypeParam{42});
    }

    { // leading + sign (or do we want this to succeed?)
        std::vector<char> const str{'+', '1', '3'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(res.ptr, &str[0]);
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
        EXPECT_EQ(value, TypeParam{42});
    }

    { // hexadecimal x prefix is not recognized
        std::vector<char> const str{'x', '3', 'f'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value, 16);

        EXPECT_EQ(res.ptr, &str[0]);
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
        EXPECT_EQ(value, TypeParam{42});
    }
}

TYPED_TEST(integral_from_char_test, binary_number)
{
    TypeParam value{42};
    std::vector<char> const str{'1', '1', '0', '1'};

    auto res = std::from_chars(&str[0], &str[0] + str.size(), value, 2);

    EXPECT_EQ(res.ptr, &str[0] + str.size());
    EXPECT_EQ(res.ec, std::errc{});
    EXPECT_EQ(value, TypeParam{13});
}

TYPED_TEST(integral_from_char_test, hexadicimal_number)
{
    TypeParam value{42};
    {
        // NOTE: According to cppreference (03.10.2018)
        // https://en.cppreference.com/w/cpp/utility/from_chars
        // uppercase letter should not be recognised so this test might need to
        // be updated when <charconv> header is included into the stl.
        std::vector<char> const str{'3', 'F'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value, 16);

        EXPECT_EQ(res.ptr, &str[0] + str.size());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{63});
    }

    value = 42; // reset
    {
        std::vector<char> const str{'3', 'f'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value, 16);

        EXPECT_EQ(res.ptr, &str[0] + str.size());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{63});
    }
}

// =============================================================================
// std::to_chars for integral types
// =============================================================================

TYPED_TEST(integral_from_char_test, to_chars)
{
    uint8_t max_num_digits = static_cast<uint8_t>(std::log10(std::numeric_limits<TypeParam>::max())) + 1;
    std::array<char, 21> buffer{}; // num_digits of uint64_t is 20 and we want to do buffer[num_digits]

    TypeParam val{0};
    for (uint8_t num_digits = 1; num_digits <= max_num_digits; ++num_digits)
    {
        // 1, 12, 123, 1234, 12345, ....
        val *= 10;
        val += num_digits % 10;
        auto res = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);

        EXPECT_EQ(res.ptr, &buffer[num_digits]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ((std::string_view{buffer.data(), num_digits}), std::string_view{std::to_string(val)});
    }
}

TYPED_TEST(integral_from_char_test, to_chars_small_value)
{
    TypeParam val{120};
    std::array<char, 10> buffer{};

    auto res = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);

    EXPECT_EQ(res.ptr, &buffer[3]);
    EXPECT_EQ(res.ec, std::errc{});
    EXPECT_EQ((std::string_view{buffer.data(), 3}), std::string_view{"120"});
}

TYPED_TEST(integral_from_char_test, to_chars_error)
{
    TypeParam val{120};
    std::array<char, 1> buffer{};

    auto res = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);

    EXPECT_EQ(res.ptr, buffer.data() + buffer.size());
    EXPECT_EQ(res.ec, std::errc::value_too_large);
}

// https://github.com/seqan/seqan3/issues/1595
TEST(to_chars_test, issue_1595)
{
    uint64_t val{123'456'789};
    std::array<char, 100> buffer{};

    auto res = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);

    EXPECT_EQ(res.ptr, &buffer[9]);
    EXPECT_EQ(res.ec, std::errc{});
    EXPECT_EQ((std::string_view{buffer.data(), 9}), std::string_view{"123456789"});
}
