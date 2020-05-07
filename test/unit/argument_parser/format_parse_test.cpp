// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/argument_parser/argument_parser.hpp>

TEST(parse_type_test, add_option_short_id)
{
    std::string option_value;

    const char * argv[] = {"./argument_parser_test", "-s", "option_string"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_section("My options");       // no-op for code coverage
    parser.add_subsection("My suboptions"); // no-op for code coverage
    parser.add_line("line");                // no-op for code coverage
    parser.add_list_item("list", "item");   // no-op for code coverage
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    // add with no space
    const char * argv2[] = {"./argument_parser_test", "-Soption_string"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, false};
    parser2.add_option(option_value, 'S', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    // ad with = sign
    const char * argv3[] = {"./argument_parser_test", "-s=option_string"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    parser3.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");
}

TEST(parse_type_test, add_option_long_id)
{
    std::string option_value;

    const char * argv[] = {"./argument_parser_test", "--string-option", "option_string"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    const char * argv3[] = {"./argument_parser_test", "--string-option=option_string"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    parser3.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");
}

TEST(parse_type_test, add_flag_short_id_single)
{
    bool option_value1{false};
    bool option_value2{true};

    const char * argv[] = {"./argument_parser_test", "-t"};
    seqan3::argument_parser parser{"test_parser", 2, argv, false};
    parser.add_flag(option_value1, 't', "true-flag", "this is a flag.");
    parser.add_flag(option_value2, 'f', "false-flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value1, true);
    EXPECT_EQ(option_value2, false);
}

TEST(parse_type_test, add_flag_short_id_multiple)
{
    bool option_value1{false};
    bool option_value2{true};
    bool option_value3{false};
    bool option_value4{false};

    const char * argv[] = {"./argument_parser_test", "-tab"};
    seqan3::argument_parser parser{"test_parser", 2, argv, false};
    parser.add_flag(option_value1, 't', "true-flag", "this is a flag.");
    parser.add_flag(option_value2, 'f', "false-flag", "this is a flag.");
    parser.add_flag(option_value3, 'a', "additional-flag", "this is a flag.");
    parser.add_flag(option_value4, 'b', "another-flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value1, true);
    EXPECT_EQ(option_value2, false);
    EXPECT_EQ(option_value3, true);
    EXPECT_EQ(option_value4, true);
}

TEST(parse_type_test, add_flag_long_id)
{
    bool option_value1{false};
    bool option_value2{true};

    const char * argv[] = {"./argument_parser_test", "--true-flag"};
    seqan3::argument_parser parser{"test_parser", 2, argv, false};
    parser.add_flag(option_value1, 't', "true-flag", "this is a flag.");
    parser.add_flag(option_value2, 'f', "false-flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value1, true);
    EXPECT_EQ(option_value2, false);
}

TEST(parse_type_test, add_positional_option)
{
    std::string positional_value;

    const char * argv[] = {"./argument_parser_test", "positional_string"};
    seqan3::argument_parser parser{"test_parser", 2, argv, false};
    parser.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "positional_string");
}

TEST(parse_type_test, independent_add_order)
{
    // testing same command line input different add_* order

    std::string positional_value;
    bool flag_value;
    int option_value;

    // Order 1: option, flag, positional
    const char * argv[] = {"./argument_parser_test", "-i", "2", "-b", "arg"};
    seqan3::argument_parser parser{"test_parser", 5, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: flag, option, positional
    seqan3::argument_parser parser2{"test_parser", 5, argv, false};
    parser2.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser2.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: option, positional, flag
    seqan3::argument_parser parser3{"test_parser", 5, argv, false};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser3.add_positional_option(positional_value, "this is a string positional.");
    parser3.add_flag(flag_value, 'b', "flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: flag, positional, option
    seqan3::argument_parser parser4{"test_parser", 5, argv, false};
    parser4.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser4.add_positional_option(positional_value, "this is a string positional.");
    parser4.add_option(option_value, 'i', "int-option", "this is a int option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: positional, flag, option
    seqan3::argument_parser parser5{"test_parser", 5, argv, false};
    parser5.add_positional_option(positional_value, "this is a string positional.");
    parser5.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser5.add_option(option_value, 'i', "int-option", "this is a int option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: positional, option, flag
    seqan3::argument_parser parser6{"test_parser", 5, argv, false};
    parser6.add_positional_option(positional_value, "this is a string positional.");
    parser6.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser6.add_flag(flag_value, 'b', "flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser6.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);
}

TEST(parse_type_test, independent_cmd_order)
{
    // testing different command line order

    std::string positional_value;
    bool flag_value;
    int option_value;

    // Order 1: option, flag, positional (POSIX conform)
    const char * argv[] = {"./argument_parser_test", "-i", "2", "-b", "arg"};
    seqan3::argument_parser parser{"test_parser", 5, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: flag, option, positional (POSIX conform)
    const char * argv2[] = {"./argument_parser_test", "-b", "-i", "2", "arg"};
    seqan3::argument_parser parser2{"test_parser", 5, argv2, false};
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser2.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser2.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: option, positional, flag
    const char * argv3[] = {"./argument_parser_test", "-i", "2", "arg", "-b"};
    seqan3::argument_parser parser3{"test_parser", 5, argv3, false};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser3.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser3.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: flag, positional, option
    const char * argv4[] = {"./argument_parser_test", "-b", "arg", "-i", "2"};
    seqan3::argument_parser parser4{"test_parser", 5, argv4, false};
    parser4.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser4.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser4.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: positional, flag, option
    const char * argv5[] = {"./argument_parser_test", "arg", "-b", "-i", "2"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, false};
    parser5.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser5.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser5.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    // Order 1: positional, option, flag
    const char * argv6[] = {"./argument_parser_test", "arg", "-i", "2", "-b"};
    seqan3::argument_parser parser6{"test_parser", 5, argv6, false};
    parser6.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser6.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser6.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser6.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);
}

TEST(parse_test, double_dash_separation_success)
{
    std::string option_value;

    // string option with dash
    const char * argv[] = {"./argument_parser_test", "--", "-strange"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_positional_option(option_value, "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "-strange");

    // negative integer option
    int option_value_int;
    const char * argv2[] = {"./argument_parser_test", "--", "-120"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_positional_option(option_value_int, "this is a int option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value_int, -120);
}

TEST(parse_test, special_characters_as_value_success)
{
    std::string option_value;

    // weird option value. But since r/regex option is parsed and with it's
    // value should work correct
    const char * argv[] = {"./argument_parser_test", "--regex", "-i=/45*&//--"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'r', "regex", "strange option value.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "-i=/45*&//--");
}

TEST(parse_test, empty_value_error)
{
    int option_value;

    // short option
    const char * argv[] = {"./argument_parser_test", "-i"};
    seqan3::argument_parser parser{"test_parser", 2, argv, false};
    parser.add_option(option_value, 'i', "long", "this is a int option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // long option
    const char * argv2[] = {"./argument_parser_test", "--long"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, false};
    parser2.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);

    // short option
    const char * argv3[] = {"./argument_parser_test", "-i="};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    parser3.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser3.parse(), seqan3::argument_parser_error);

    // short option
    const char * argv4[] = {"./argument_parser_test", "--long="};
    seqan3::argument_parser parser4{"test_parser", 2, argv4, false};
    parser4.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser4.parse(), seqan3::argument_parser_error);
}

TEST(parse_type_test, parse_success_bool_option)
{
    bool option_value;
    bool positional_value;

    // numbers 0 and 1
    {
        const char * argv[] = {"./argument_parser_test", "-b", "1", "0"};
        seqan3::argument_parser parser{"test_parser", 4, argv, false};
        parser.add_option(option_value, 'b', "bool-option", "this is a bool option.");
        parser.add_positional_option(positional_value, "this is a bool positional.");

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, true);
        EXPECT_EQ(positional_value, false);
    }

    // true and false
    {
        const char * argv[] = {"./argument_parser_test", "-b", "true", "false"};
        seqan3::argument_parser parser{"test_parser", 4, argv, false};
        parser.add_option(option_value, 'b', "bool-option", "this is a bool option.");
        parser.add_positional_option(positional_value, "this is a bool positional.");

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, true);
        EXPECT_EQ(positional_value, false);
    }
}

TEST(parse_type_test, parse_success_int_option)
{
    int option_value;
    size_t positional_value;

    const char * argv[] = {"./argument_parser_test", "-i", "-2", "278"};
    seqan3::argument_parser parser{"test_parser", 4, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser.add_positional_option(positional_value, "this is a int positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());

    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -2);
    EXPECT_EQ(positional_value, 278u);
}

TEST(parse_type_test, parse_success_double_option)
{
    double option_value;
    double positional_value;

    const char * argv[] = {"./argument_parser_test", "-d", "12.457", "0.123"};
    seqan3::argument_parser parser{"test_parser", 4, argv, false};
    parser.add_option(option_value, 'd', "double-option", "this is a double option.");
    parser.add_positional_option(positional_value, "this is a double positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());

    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_FLOAT_EQ(option_value, 12.457);
    EXPECT_FLOAT_EQ(positional_value, 0.123);

    // double expression with e
    const char * argv2[] = {"./argument_parser_test", "-d", "6.0221418e23"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'd', "double-option", "this is a double option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());

    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_FLOAT_EQ(option_value, 6.0221418e23);
    EXPECT_FLOAT_EQ(positional_value, 0.123);

}

TEST(parse_type_test, parse_error_bool_option)
{
    bool option_value;

    // fail on character input
    const char * argv[] = {"./argument_parser_test", "-b", "a"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'b', "bool-option", "this is a bool option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // fail on number input expect 0 and 1
    const char * argv2[] = {"./argument_parser_test", "-b", "124"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'b', "bool-option", "this is a bool option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);
}

TEST(parse_type_test, parse_error_int_option)
{
    int option_value;

    // fail on character
    const char * argv[] = {"./argument_parser_test", "-i", "abc"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // fail on number followed by character
    const char * argv2[] = {"./argument_parser_test", "-i", "2abc"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);

    // fail on double
    const char * argv3[] = {"./argument_parser_test", "-i", "3.12"};
    seqan3::argument_parser parser3{"test_parser", 3, argv3, false};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser3.parse(), seqan3::argument_parser_error);

    // fail on negative number for unsigned
    unsigned option_value_u;
    const char * argv4[] = {"./argument_parser_test", "-i", "-1"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, false};
    parser4.add_option(option_value_u, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser4.parse(), seqan3::argument_parser_error);

    // fail on overflow
    int8_t option_value_int8;
    const char * argv5[] = {"./argument_parser_test", "-i", "129"};
    seqan3::argument_parser parser5{"test_parser", 3, argv5, false};
    parser5.add_option(option_value_int8, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser5.parse(), seqan3::argument_parser_error);

    uint8_t option_value_uint8;
    const char * argv6[] = {"./argument_parser_test", "-i", "267"};
    seqan3::argument_parser parser6{"test_parser", 3, argv6, false};
    parser6.add_option(option_value_uint8, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser6.parse(), seqan3::argument_parser_error);
}

TEST(parse_type_test, parse_error_double_option)
{
    double option_value;

    // fail on character
    const char * argv[] = {"./argument_parser_test", "-i", "abc"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'd', "double-option", "this is a double option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // fail on number followed by character
    const char * argv2[] = {"./argument_parser_test", "-d", "12.457a"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'd', "double-option", "this is a double option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);
}

namespace foo
{
enum class bar
{
    one,
    two,
    three
};

auto enumeration_names(bar)
{
    return std::unordered_map<std::string_view, bar>{{"one", bar::one}, {"two", bar::two}, {"three", bar::three}};
}
} // namespace foo

namespace Other
{
enum class bar
{
    one,
    two
};
} // namespace Other

namespace seqan3::custom
{
template <>
struct argument_parsing<Other::bar>
{
    static inline std::unordered_map<std::string_view, Other::bar> const enumeration_names
    {
        {"one", Other::bar::one}, {"two", Other::bar::two}
    };
};
} // namespace seqan3::custom

TEST(parse_type_test, parse_success_enum_option)
{
    {
        foo::bar option_value{};

        const char * argv[] = {"./argument_parser_test", "-e", "two"};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE(option_value == foo::bar::two);
    }

    {
        Other::bar option_value{};

        const char * argv[] = {"./argument_parser_test", "-e", "two"};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE(option_value == Other::bar::two);
    }
}

TEST(parse_type_test, parse_error_enum_option)
{
    foo::bar option_value{};

    const char * argv[] = {"./argument_parser_test", "-e", "four"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);
}

TEST(parse_test, too_many_arguments_error)
{
    int option_value;

    const char * argv[] = {"./argument_parser_test", "5", "15"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    parser.add_positional_option(option_value, "this is an int option.");

    EXPECT_THROW(parser.parse(), seqan3::too_many_arguments);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv2[] = {"./argument_parser_test", "2", "--", "-i"};
    seqan3::argument_parser parser2{"test_parser", 4, argv2, false};
    parser2.add_positional_option(option_value, "normal int positional argument.");
    parser2.add_option(option_value, 'i', "int-option", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::too_many_arguments);
}

TEST(parse_test, too_few_arguments_error)
{
    int option_value;

    const char * argv[] = {"./argument_parser_test", "15"};
    seqan3::argument_parser parser{"test_parser", 2, argv, false};
    parser.add_positional_option(option_value, "this is an int option.");
    parser.add_positional_option(option_value, "this is another option.");

    EXPECT_THROW(parser.parse(), seqan3::too_few_arguments);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv2[] = {"./argument_parser_test", "-i", "2"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_positional_option(option_value, "normal int positional argument.");
    parser2.add_option(option_value, 'i', "int-option", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::too_few_arguments);
}

TEST(parse_test, unknown_option_error)
{
    // unknown short option
    const char * argv[] = {"./argument_parser_test", "-i", "15"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};

    EXPECT_THROW(parser.parse(), seqan3::unknown_option);

    // unknown long option
    const char * argv2[] = {"./argument_parser_test", "--arg", "8"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};

    EXPECT_THROW(parser2.parse(), seqan3::unknown_option);

    // unknown short flag
    const char * argv3[] = {"./argument_parser_test", "-a"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};

    EXPECT_THROW(parser3.parse(), seqan3::unknown_option);

    // unknown long flag
    const char * argv4[] = {"./argument_parser_test", "--arg"};
    seqan3::argument_parser parser4{"test_parser", 2, argv4, false};

    EXPECT_THROW(parser4.parse(), seqan3::unknown_option);

    // negative numbers are seen as options
    const char * argv5[] = {"./argument_parser_test", "-5"};
    seqan3::argument_parser parser5{"test_parser", 2, argv5, false};

    EXPECT_THROW(parser5.parse(), seqan3::unknown_option);

    // unknown short option in more complex command line
    int option_value_i;
    std::string option_value_a;
    std::string positional_option;
    const char * argv6[] = {"./argument_parser_test", "-i", "129", "arg1", "-b", "bcd", "-a", "abc"};
    seqan3::argument_parser parser6{"test_parser", 8, argv6, false};
    parser6.add_option(option_value_i, 'i', "int-option", "this is a int option.");
    parser6.add_option(option_value_a, 'a', "string-option", "this is a string option.");
    parser6.add_positional_option(positional_option, "normal int positional argument.");

    EXPECT_THROW(parser6.parse(), seqan3::unknown_option);
}

TEST(parse_test, option_declared_multiple_times_error)
{
    int option_value;

    // short option
    const char * argv[] = {"./argument_parser_test", "-i", "15", "-i", "3"};
    seqan3::argument_parser parser{"test_parser", 5, argv, false};
    parser.add_option(option_value, 'i', "long", "this is a int option.");

    EXPECT_THROW(parser.parse(), seqan3::option_declared_multiple_times);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv2[] = {"./argument_parser_test", "--long", "5", "--long", "6"};
    seqan3::argument_parser parser2{"test_parser", 5, argv2, false};
    parser2.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::option_declared_multiple_times);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv3[] = {"./argument_parser_test", "-i", "5", "--long", "6"};
    seqan3::argument_parser parser3{"test_parser", 5, argv3, false};
    parser3.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser3.parse(), seqan3::option_declared_multiple_times);
}

TEST(parse_test, required_option_missing)
{
    int option_value;

    // option is required
    const char * argv[] = {"./argument_parser_test", "5", "-i", "15"};
    seqan3::argument_parser parser{"test_parser", 4, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is an int option.");
    parser.add_option(option_value, 'a', "req-option", "I am required.", seqan3::option_spec::REQUIRED);
    parser.add_positional_option(option_value, "this is an int option.");

    EXPECT_THROW(parser.parse(), seqan3::required_option_missing);
}

TEST(parse_test, argv_const_combinations)
{
    bool flag_value{false};

    char arg1[]{"./argument_parser"};
    char arg2[]{"-f"};
    char * argv[] = {arg1, arg2};

    // all const*
    char const * const * const argv_all_const{argv};
    seqan3::argument_parser parser{"test_parser", 2, argv_all_const, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // none const
    flag_value = false;
    parser = seqan3::argument_parser{"test_parser", 2, argv, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 1
    flag_value = false;
    char const * argv_const1[] = {"./argument_parser_test", "-f"};
    parser = seqan3::argument_parser{"test_parser", 2, argv_const1, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 2
    flag_value = false;
    char * const argv_const2[] = {arg1, arg2};
    parser = seqan3::argument_parser{"test_parser", 2, argv_const2, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 12
    flag_value = false;
    char const * const argv_const12[] = {arg1, arg2};
    parser = seqan3::argument_parser{"test_parser", 2, argv_const12, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);
}

TEST(parse_test, multiple_empty_options)
{
    int option_value;

    {
        const char * argv[]{"./empty_long", "-s=1"};
        seqan3::argument_parser parser{"empty_long", 2, argv, false};
        parser.add_option(option_value, 'i', "", "no long");
        parser.add_option(option_value, 's', "", "no long");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(1, option_value);
    }

    {
        const char * argv[]{"./empty_long", "-s=1", "--unknown"};
        seqan3::argument_parser parser{"empty_long", 3, argv, false};
        parser.add_option(option_value, 'i', "", "no long");
        parser.add_option(option_value, 's', "", "no long");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    {
        const char * argv[]{"./empty_short", "--long=2"};
        seqan3::argument_parser parser{"empty_short", 2, argv, false};
        parser.add_option(option_value, '\0', "longi", "no short");
        parser.add_option(option_value, '\0', "long", "no short");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(2, option_value);
    }
}

TEST(parse_test, version_check_option_error)
{
    {   // version-check must be followed by a value
        const char * argv[] = {"./argument_parser_test", "--version-check"};
        EXPECT_THROW((seqan3::argument_parser{"test_parser", 2, argv}), seqan3::argument_parser_error);
    }

    {   // version-check value must be 0 or 1
        const char * argv[] = {"./argument_parser_test", "--version-check", "foo"};
        EXPECT_THROW((seqan3::argument_parser{"test_parser", 3, argv}), seqan3::argument_parser_error);
    }
}

TEST(parse_test, subcommand_argument_parser_success)
{
    bool flag_value{};
    std::string option_value{};

    // parsing
    {
        const char * argv[]{"./top_level", "-f", "sub1", "foo"};
        seqan3::argument_parser top_level_parser{"top_level", 4, argv, false, {"sub1", "sub2"}};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        EXPECT_NO_THROW(top_level_parser.parse());
        EXPECT_EQ(true, flag_value);

        seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();

        EXPECT_EQ(sub_parser.info.app_name, "top_level-sub1");

        sub_parser.add_positional_option(option_value, "foo bar");

        EXPECT_NO_THROW(sub_parser.parse());
        EXPECT_EQ("foo", option_value);
    }

    // top-level help page
    {
        const char * argv[]{"./top_level", "-h", "-f", "sub1", "foo"};
        seqan3::argument_parser top_level_parser{"top_level", 5, argv, false, {"sub1", "sub2"}};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        testing::internal::CaptureStdout();
        EXPECT_EXIT(top_level_parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        EXPECT_FALSE(std::string{testing::internal::GetCapturedStdout()}.empty());
    }

    // sub-parser help page
    {
        const char * argv[]{"./top_level", "-f", "sub1", "-h"};
        seqan3::argument_parser top_level_parser{"top_level", 4, argv, false, {"sub1", "sub2"}};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        EXPECT_NO_THROW(top_level_parser.parse());
        EXPECT_EQ(true, flag_value);

        seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();

        EXPECT_EQ(sub_parser.info.app_name, "top_level-sub1");

        sub_parser.add_positional_option(option_value, "foo bar");

        testing::internal::CaptureStdout();
        EXPECT_EXIT(sub_parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        EXPECT_FALSE(std::string{testing::internal::GetCapturedStdout()}.empty());
    }

    // incorrect sub command
    {
        const char * argv[]{"./top_level", "-f", "2", "subiddysub", "foo"};
        EXPECT_THROW((seqan3::argument_parser{"top_level", 5, argv, false, {"sub1", "sub2"}}),
                     seqan3::argument_parser_error);
    }
}

TEST(parse_test, issue1544)
{
    {   // wrong separation of long value:
        std::string option_value;
        const char * argv[] = {"./argument_parser_test", "--foohallo"};
        seqan3::argument_parser parser{"test_parser", 2, argv, false};
        parser.add_option(option_value, 'f', "foo", "this is a string option.");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    {   // unknown option (`--foo-bar`) that has a prefix of a known option (`--foo`)
        std::string option_value;
        const char * argv[] = {"./argument_parser_test", "--foo", "hallo", "--foo-bar", "ballo"};
        seqan3::argument_parser parser{"test_parser", 5, argv, false};
        parser.add_option(option_value, 'f', "foo", "this is a string option.");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    {   // known option (`--foo-bar`) that has a prefix of a unknown option (`--foo`)
        std::string option_value;
        const char * argv[] = {"./argument_parser_test", "--foo", "hallo", "--foo-bar", "ballo"};
        seqan3::argument_parser parser{"test_parser", 5, argv, false};
        parser.add_option(option_value, 'f', "foo-bar", "this is a string option.");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    {   // known option (`--foo`) is a prefix of another known option (`--foo-bar`)
        std::string foo_option_value;
        std::string foobar_option_value;
        const char * argv[] = {"./argument_parser_test", "--foo", "hallo", "--foo-bar", "ballo"};
        seqan3::argument_parser parser{"test_parser", 5, argv, false};
        parser.add_option(foo_option_value, 'f', "foo", "this is a prefix of foobar.");
        parser.add_option(foobar_option_value, 'b', "foo-bar", "this has prefix foo.");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(foo_option_value, "hallo");
        EXPECT_EQ(foobar_option_value, "ballo");
    }
}
