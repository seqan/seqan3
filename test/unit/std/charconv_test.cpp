// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <limits>

#include <seqan3/std/charconv>
#include <seqan3/std/concepts>

// =============================================================================
// std::from_chars for integral types
// =============================================================================

template <typename T>
class integral_from_char_test: public ::testing::Test { };

using integral_types = ::testing::Types<int8_t, uint8_t, int16_t, uint16_t,
                                        int32_t, uint32_t, int64_t, uint64_t>;

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
    std::vector<char> const str{'1', '2', '3', '0', '0', '0', '0', '0', '0', '0',
                                '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0'};

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
    { // trailing char
        std::vector<char> const str{'1', '2', 'a'};

        auto res = std::from_chars(&str[0], &str[0] + str.size(), value);

        EXPECT_EQ(res.ptr, &str[2]);
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(value, TypeParam{12});
    }

    value = 42; // reset
    { // float
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

// =============================================================================
// std::from_chars for float, double and long double
// =============================================================================

template <typename T>
class from_char_real_test: public ::testing::Test { };

using real_types = ::testing::Types<float, double, long double>;

TYPED_TEST_SUITE(from_char_real_test, real_types, );

TYPED_TEST(from_char_real_test, real_numbers)
{
    std::setlocale(LC_NUMERIC, "C");
    {
        TypeParam val{};
        std::string str = "1234";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{1234});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "1.2e3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{1200});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "1.2e-3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{0.0012});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "1.e2";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{100});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "1.";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{1});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = ".2e3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{200});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "2e3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{2000});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "2";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{2});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "4em";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{4});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + 1);
    }

    {
        TypeParam val{};
        std::string str = "-1.2e3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{-1200});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{42};
        std::string str = "-.3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{-0.3});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{42};
        std::string str = "1.2e";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{1.2});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + 3);
    }

    {
        TypeParam val{42};
        std::string str = "0.0";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{0});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    // Read only until a certain position
    {
        TypeParam val{42};
        std::string str = "3.194357";
        auto res = std::from_chars(&str[0], &str[0] + 4, val);
        EXPECT_FLOAT_EQ(val, TypeParam{3.19});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + 4);
    }

    // Partial Parsing
    {
        TypeParam val{42};
        std::string str = "3.19abc";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{3.19});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + 4);
    }
}

TYPED_TEST(from_char_real_test, infinity_value)
{
    {
        TypeParam val{};
        std::string str = "inf";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_EQ(val, std::numeric_limits<TypeParam>::infinity());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "infinity";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_EQ(val, std::numeric_limits<TypeParam>::infinity());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "INF";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_EQ(val, std::numeric_limits<TypeParam>::infinity());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "INFINITY";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_EQ(val, std::numeric_limits<TypeParam>::infinity());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

}

TYPED_TEST(from_char_real_test, nan_value)
{
    // Note:
    // According to the IEEE standard, NaN values have the odd property that
    // comparisons involving them are always false. That is, for a float f,
    // f != f will be true only if f is NaN.
    {
        TypeParam val{};
        std::string str = "nan";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_TRUE(std::isnan(val));
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "NAN";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_TRUE(std::isnan(val));
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "nan(abc)";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_TRUE(std::isnan(val));
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }

    {
        TypeParam val{};
        std::string str = "NAN(abc)";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_TRUE(std::isnan(val));
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, &str[0] + str.size());
    }
}

TYPED_TEST(from_char_real_test, non_valid_strings)
{
    {
        TypeParam val{42};
        std::string str = "e3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{42});
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
    }

    {
        TypeParam val{42};
        std::string str = "+1.2e3";
        auto res = std::from_chars(&str[0], &str[0] + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{42});
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
    }
}

// =============================================================================
// std::to_chars for float, double and long double
// =============================================================================

TYPED_TEST(from_char_real_test, to_chars)
{
    TypeParam val{120.3};
    std::array<char, 10> buffer{};

    auto res = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);

    EXPECT_EQ(res.ptr, &buffer[5]);
    EXPECT_EQ(res.ec, std::errc{});
    EXPECT_EQ((std::string_view{buffer.data(), 5}), std::string_view{"120.3"});
}
