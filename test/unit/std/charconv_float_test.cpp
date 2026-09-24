// SPDX-FileCopyrightText: 2006-2026 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2026 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

// make sure that including the std header does not produce any errors
// see https://github.com/seqan/seqan3/issues/2352
#include <seqan3/std/charconv>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <iostream>
#include <limits>
#include <string_view>

// =============================================================================
// std::from_chars for float, double and long double
// =============================================================================

template <typename T>
class from_char_real_test : public ::testing::Test
{};

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
// from_chars with chars_format and partial ranges
// =============================================================================

// "1234567890" repeated ten times and "8999" at the end.
static constexpr std::string_view long_number{"123456789012345678901234567890123456789012345678901234567890"
                                              "12345678901234567890123456789012345678999"};
static_assert(long_number.size() == 101u);

TYPED_TEST(from_char_real_test, general_partial_range)
{
    // The range [first, last) excludes the last character.
    {
        TypeParam val{42};
        auto res = std::from_chars(long_number.data(),
                                   long_number.data() + long_number.size() - 1u,
                                   val,
                                   std::chars_format::general);
        EXPECT_EQ(res.ptr, long_number.data() + long_number.size() - 1u);
        if constexpr (std::is_same_v<TypeParam, float>)
        {
            EXPECT_EQ(res.ec, std::errc::result_out_of_range); // libc++ sets val to inf.
        }
        else
        {
            EXPECT_EQ(res.ec, std::errc{});
            EXPECT_NEAR(val, 1.2345678901234567890e99, 1e84);
        }
    }

    {
        TypeParam val{42};
        std::string_view str{"15 foo"};
        auto res = std::from_chars(str.data(), str.data() + str.size() - 1u, val, std::chars_format::general);
        EXPECT_FLOAT_EQ(val, TypeParam{15});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + 2);
    }

    {
        TypeParam val{42};
        std::string_view str{"5000000000"};
        auto res = std::from_chars(str.data(), str.data() + str.size() - 1u, val, std::chars_format::general);
        EXPECT_FLOAT_EQ(val, TypeParam{500'000'000});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size() - 1u);
    }

    // Longer than the internal stack buffer.
    {
        TypeParam val{42};
        std::string str = std::string(150u, '0') + "12.5e1";
        auto res = std::from_chars(str.data(), str.data() + str.size() - 1u, val, std::chars_format::general);
        EXPECT_FLOAT_EQ(val, TypeParam{12.5});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size() - 2u);
    }

    for (std::string_view const str : {"bar", " 42", "+42"})
    {
        TypeParam val{42};
        auto res = std::from_chars(str.data(), str.data() + str.size() - 1u, val, std::chars_format::general);
        EXPECT_FLOAT_EQ(val, TypeParam{42}) << str;
        EXPECT_EQ(res.ec, std::errc::invalid_argument) << str;
        EXPECT_EQ(res.ptr, str.data()) << str;
    }

    // The 0x prefix is not permitted: "0x123" is parsed as the value "0" with unparsed remainder "x123".
    {
        TypeParam val{42};
        std::string_view str{"0x123"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{0});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + 1);
    }

    // A stale errno does not affect the result.
    {
        TypeParam val{42};
        std::string_view str{"1.5"};
        errno = ERANGE;
        auto res = std::from_chars(str.data(), str.data() + str.size(), val);
        EXPECT_FLOAT_EQ(val, TypeParam{1.5});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size());
    }
}

TYPED_TEST(from_char_real_test, hex)
{
    {
        TypeParam val{42};
        std::string_view str{""};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{42});
        EXPECT_EQ(res.ec, std::errc::invalid_argument);
        EXPECT_EQ(res.ptr, str.data());
    }

    {
        TypeParam val{42};
        auto res =
            std::from_chars(long_number.data(), long_number.data() + long_number.size(), val, std::chars_format::hex);
        EXPECT_EQ(res.ptr, long_number.data() + long_number.size());
        if constexpr (std::is_same_v<TypeParam, float>)
        {
            EXPECT_EQ(res.ec, std::errc::result_out_of_range); // libc++ sets val to inf.
        }
        else
        {
            EXPECT_EQ(res.ec, std::errc{});
            EXPECT_NEAR(val, 0x1.2345678901234567890p400, 1e105);
        }
    }

    {
        TypeParam val{42};
        std::string_view str{"15 foo"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{0x15});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + 2);
    }

    {
        TypeParam val{42};
        std::string_view str{"bar"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{0xba});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + 2);
    }

    {
        TypeParam val{42};
        std::string_view str{"5000000000"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{0x50'00'00'00'00});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size());
    }

    {
        TypeParam val{42};
        std::string_view str{"-1a"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{-0x1a});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size());
    }

    // The 0x prefix is not permitted: "0x123" is parsed as the value "0" with unparsed remainder "x123".
    {
        TypeParam val{42};
        std::string_view str{"0x123"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{0});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + 1);
    }

    {
        TypeParam val{};
        std::string_view str{"-inf"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_EQ(val, -std::numeric_limits<TypeParam>::infinity());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size());
    }

    for (std::string_view const str : {" 42", "+42", "xyz"})
    {
        TypeParam val{42};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::hex);
        EXPECT_FLOAT_EQ(val, TypeParam{42}) << str;
        EXPECT_EQ(res.ec, std::errc::invalid_argument) << str;
        EXPECT_EQ(res.ptr, str.data()) << str;
    }
}

TYPED_TEST(from_char_real_test, fixed_and_scientific)
{
    // fixed: The exponent is not part of the pattern.
    {
        TypeParam val{42};
        std::string_view str{"1.5e3"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::fixed);
        EXPECT_FLOAT_EQ(val, TypeParam{1.5});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + 3);
    }

    // scientific: The exponent is required.
    {
        TypeParam val{42};
        std::string_view str{"1.5e3"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::scientific);
        EXPECT_FLOAT_EQ(val, TypeParam{1500});
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size());
    }

    for (std::string_view const str : {"1.5", "1.5e", "1.5e+"})
    {
        TypeParam val{42};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::scientific);
        EXPECT_FLOAT_EQ(val, TypeParam{42}) << str;
        EXPECT_EQ(res.ec, std::errc::invalid_argument) << str;
        EXPECT_EQ(res.ptr, str.data()) << str;
    }

    {
        TypeParam val{};
        std::string_view str{"inf"};
        auto res = std::from_chars(str.data(), str.data() + str.size(), val, std::chars_format::scientific);
        EXPECT_EQ(val, std::numeric_limits<TypeParam>::infinity());
        EXPECT_EQ(res.ec, std::errc{});
        EXPECT_EQ(res.ptr, str.data() + str.size());
    }
}

// =============================================================================
// std::to_chars for float, double and long double
// =============================================================================

TYPED_TEST(from_char_real_test, to_chars)
{
    // We use a power of two (i.e. 2^(-2)) for the fractional part to have a stable floating point number across
    // different floating point types (e.g. `float`, `double`, `long double`).
    // Other values, lets say 120.3, could have different string representations like `120.3`, `120.30000...01`, or
    // `120.2999...9716` depending on the actual implementation.
    TypeParam val{120.25};
    std::array<char, 10> buffer{};

    auto res = std::to_chars(buffer.data(), buffer.data() + buffer.size(), val);
    size_t used_buffer_size = res.ptr - buffer.data();

    EXPECT_EQ(used_buffer_size, 6u);
    EXPECT_EQ(res.ec, std::errc{});
    EXPECT_EQ((std::string_view{buffer.data(), used_buffer_size}), std::string_view{"120.25"});
}
