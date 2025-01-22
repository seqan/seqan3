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
