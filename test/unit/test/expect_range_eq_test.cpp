// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h> // provides test utility to test google test itself

#include <seqan3/std/span>
#include <string_view>

#include <seqan3/test/expect_range_eq.hpp>

TEST(expect_range_eq, braces_with_many_commas)
{
    std::vector expect{0, 1, 2};
    EXPECT_RANGE_EQ(expect, std::vector({0, 1, 2}));

    // Note: the macro can't handle the following expression,
    //     EXPECT_RANGE_EQ(expect, std::vector{0, 1, 2});
    // because it confuses it with a call with 4 arguments:
    //     EXPECT_RANGE_EQ("expect", "std::vector{0", "1", "2"});

    // workaround:
    EXPECT_RANGE_EQ(expect, (std::vector{0, 1, 2}));
}

TEST(string_view, range_eq_pass)
{
    std::vector<char> expect{'H', 'e', 'l', 'l', 'o'};
    std::string_view result{"Hello"};

    auto && expect_result = seqan3::test::expect_range_eq{}("expect", "result", expect, result);
    EXPECT_TRUE(expect_result);
    EXPECT_RANGE_EQ(expect, result);
}

TEST(string_view, range_eq_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  expect\n"
                                 "    Which is: Hel\n"
                                 "lo\n"
                                 "  result\n"
                                 "    Which is: Hello!";

    std::vector<char> expect{'H', 'e', 'l', '\n', 'l', 'o'};
    std::string_view result{"Hello!"};

    auto && expect_result = seqan3::test::expect_range_eq{}("expect", "result", expect, result);

    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(EXPECT_RANGE_EQ(expect, result), error_message);
}

TEST(span, range_eq_pass)
{
    std::vector<int> expect{0, 1, 2, 3, 4};
    std::vector<int> source{-2, -1, 0, 1, 2, 3, 4, 5, 6};
    std::span result{source.begin() + 2, 5};

    auto && expect_result = seqan3::test::expect_range_eq{}("expect", "result", expect, result);
    EXPECT_TRUE(expect_result);
    EXPECT_RANGE_EQ(expect, result);
}

TEST(span, range_eq_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  expect\n"
                                 "    Which is: [0,1,2,3,4]\n"
                                 "  result\n"
                                 "    Which is: [-1,0,1,2,3,4,5]";

    std::vector<int> expect{0, 1, 2, 3, 4};
    std::vector<int> source{-2, -1, 0, 1, 2, 3, 4, 5, 6};
    std::span result{source.begin() + 1, 7};

    auto && expect_result = seqan3::test::expect_range_eq{}("expect", "result", expect, result);

    EXPECT_FALSE(expect_result);
    EXPECT_STREQ(error_message, expect_result.message());
    EXPECT_NONFATAL_FAILURE(EXPECT_RANGE_EQ(expect, result), error_message);
}

struct input_range
{
    static constexpr int values[]{0, 1, 2, 3, 4};
    int const * current = values;
    int const * sentinel = values + sizeof(values) / sizeof(int);

    struct iterator
    {
        using difference_type = std::ptrdiff_t;
        using value_type = int;

        input_range * host;

        int const & operator*() const { return *host->current; }
        iterator & operator++() { ++host->current; return *this; }
        value_type operator++(int) { value_type x = *(*this); ++*this; return x; }
        bool operator==(iterator const &) const { return host->current == host->sentinel; }
        bool operator!=(iterator const & sentinel) const { return !(*this == sentinel);}
    };

    iterator begin() { return {this}; }
    iterator end() { return {this}; }
};

TEST(input_range, range_eq_pass)
{
    EXPECT_TRUE(std::ranges::input_range<input_range>);

    std::vector<int> expect{0, 1, 2, 3, 4};

    {
        input_range result{};
        auto && expect_result = seqan3::test::expect_range_eq{}("expect", "result", expect, result);
        EXPECT_TRUE(expect_result);
    }

    {
        input_range result{};
        EXPECT_RANGE_EQ(expect, result);
    }
}

TEST(input_range, range_eq_fail)
{
    char const * error_message = "Expected equality of these values:\n"
                                 "  expect\n"
                                 "    Which is: [0,1,2,3,4,5]\n"
                                 "  result\n"
                                 "    Which is: [0,1,2,3,4]";

    std::vector<int> expect{0, 1, 2, 3, 4, 5};

    {
        input_range result{};
        auto && expect_result = seqan3::test::expect_range_eq{}("expect", "result", expect, result);
        EXPECT_FALSE(expect_result);
        EXPECT_STREQ(error_message, expect_result.message());
    }

    {
        input_range result{};
        EXPECT_NONFATAL_FAILURE(EXPECT_RANGE_EQ(expect, result), error_message);
    }
}
