// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/argument_parser/argument_parser.hpp>

TEST(parse_type_test, add_option_short_id)
{
    std::string option_value;

    char const * argv[] = {"./argument_parser_test", "-s", "option_string"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
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
    char const * argv2[] = {"./argument_parser_test", "-Soption_string"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'S', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    // ad with = sign
    char const * argv3[] = {"./argument_parser_test", "-s=option_string"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    parser3.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");
}

TEST(parse_type_test, add_option_long_id)
{
    std::string option_value;

    char const * argv[] = {"./argument_parser_test", "--string-option", "option_string"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    char const * argv3[] = {"./argument_parser_test", "--string-option=option_string"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    parser3.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");
}

TEST(parse_type_test, add_flag_short_id_single)
{
    bool option_value1{false};
    bool option_value2{false};

    char const * argv[] = {"./argument_parser_test", "-a"};
    seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
    parser.add_flag(option_value1, 'f', "flag", "this is a flag.");
    parser.add_flag(option_value2, 'a', "another-flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value1, false);
    EXPECT_EQ(option_value2, true);
}

TEST(parse_type_test, add_flag_short_id_multiple)
{
    bool option_value1{false};
    bool option_value2{false};
    bool option_value3{false};
    bool option_value4{false};

    char const * argv[] = {"./argument_parser_test", "-acd"};
    seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
    parser.add_flag(option_value1, 'a', "flag", "this is a flag.");
    parser.add_flag(option_value2, 'b', "also-flag", "this is a flag.");
    parser.add_flag(option_value3, 'c', "additional-flag", "this is a flag.");
    parser.add_flag(option_value4, 'd', "another-flag", "this is a flag.");

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
    bool option_value2{false};

    char const * argv[] = {"./argument_parser_test", "--another-flag"};
    seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
    parser.add_flag(option_value1, 't', "flag", "this is a flag.");
    parser.add_flag(option_value2, 'f', "another-flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value1, false);
    EXPECT_EQ(option_value2, true);
}

TEST(parse_type_test, add_positional_option)
{
    std::string positional_value;

    char const * argv[] = {"./argument_parser_test", "positional_string"};
    seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
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
    bool flag_value{false};
    int option_value;

    // Order 1: option, flag, positional
    char const * argv[] = {"./argument_parser_test", "-i", "2", "-b", "arg"};
    seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false; // reinstate to default value

    // Order 1: flag, option, positional
    seqan3::argument_parser parser2{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser2.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser2.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: option, positional, flag
    seqan3::argument_parser parser3{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser3.add_positional_option(positional_value, "this is a string positional.");
    parser3.add_flag(flag_value, 'b', "flag", "this is a flag.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: flag, positional, option
    seqan3::argument_parser parser4{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser4.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser4.add_positional_option(positional_value, "this is a string positional.");
    parser4.add_option(option_value, 'i', "int-option", "this is a int option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: positional, flag, option
    seqan3::argument_parser parser5{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser5.add_positional_option(positional_value, "this is a string positional.");
    parser5.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser5.add_option(option_value, 'i', "int-option", "this is a int option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: positional, option, flag
    seqan3::argument_parser parser6{"test_parser", 5, argv, seqan3::update_notifications::off};
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
    bool flag_value{false};
    int option_value;

    // Order 1: option, flag, positional (POSIX conform)
    char const * argv[] = {"./argument_parser_test", "-i", "2", "-b", "arg"};
    seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false; // reinstate to default value

    // Order 1: flag, option, positional (POSIX conform)
    char const * argv2[] = {"./argument_parser_test", "-b", "-i", "2", "arg"};
    seqan3::argument_parser parser2{"test_parser", 5, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser2.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser2.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: option, positional, flag
    char const * argv3[] = {"./argument_parser_test", "-i", "2", "arg", "-b"};
    seqan3::argument_parser parser3{"test_parser", 5, argv3, seqan3::update_notifications::off};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser3.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser3.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: flag, positional, option
    char const * argv4[] = {"./argument_parser_test", "-b", "arg", "-i", "2"};
    seqan3::argument_parser parser4{"test_parser", 5, argv4, seqan3::update_notifications::off};
    parser4.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser4.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser4.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: positional, flag, option
    char const * argv5[] = {"./argument_parser_test", "arg", "-b", "-i", "2"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, seqan3::update_notifications::off};
    parser5.add_option(option_value, 'i', "int-option", "this is a int option.");
    parser5.add_flag(flag_value, 'b', "flag", "this is a flag.");
    parser5.add_positional_option(positional_value, "this is a string positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(positional_value, "arg");
    EXPECT_EQ(option_value, 2);
    EXPECT_EQ(flag_value, true);

    flag_value = false;

    // Order 1: positional, option, flag
    char const * argv6[] = {"./argument_parser_test", "arg", "-i", "2", "-b"};
    seqan3::argument_parser parser6{"test_parser", 5, argv6, seqan3::update_notifications::off};
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
    char const * argv[] = {"./argument_parser_test", "--", "-strange"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_positional_option(option_value, "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "-strange");

    // negative integer option
    int option_value_int;
    char const * argv2[] = {"./argument_parser_test", "--", "-120"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
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
    char const * argv[] = {"./argument_parser_test", "--regex", "-i=/45*&//--"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
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
    char const * argv[] = {"./argument_parser_test", "-i"};
    seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'i', "long", "this is a int option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // long option
    char const * argv2[] = {"./argument_parser_test", "--long"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);

    // short option
    char const * argv3[] = {"./argument_parser_test", "-i="};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    parser3.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser3.parse(), seqan3::argument_parser_error);

    // short option
    char const * argv4[] = {"./argument_parser_test", "--long="};
    seqan3::argument_parser parser4{"test_parser", 2, argv4, seqan3::update_notifications::off};
    parser4.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser4.parse(), seqan3::argument_parser_error);
}

TEST(parse_type_test, parse_success_bool_option)
{
    bool option_value;
    bool positional_value;

    // numbers 0 and 1
    {
        char const * argv[] = {"./argument_parser_test", "-b", "1", "0"};
        seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
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
        char const * argv[] = {"./argument_parser_test", "-b", "true", "false"};
        seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
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

    char const * argv[] = {"./argument_parser_test", "-i", "-2", "278"};
    seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
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

    char const * argv[] = {"./argument_parser_test", "-d", "12.457", "0.123"};
    seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'd', "double-option", "this is a double option.");
    parser.add_positional_option(positional_value, "this is a double positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());

    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_FLOAT_EQ(option_value, 12.457);
    EXPECT_FLOAT_EQ(positional_value, 0.123);

    // double expression with e
    char const * argv2[] = {"./argument_parser_test", "-d", "6.0221418e23"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
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
    char const * argv[] = {"./argument_parser_test", "-b", "a"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'b', "bool-option", "this is a bool option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // fail on number input expect 0 and 1
    char const * argv2[] = {"./argument_parser_test", "-b", "124"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'b', "bool-option", "this is a bool option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);
}

TEST(parse_type_test, parse_error_int_option)
{
    int option_value;

    // fail on character
    char const * argv[] = {"./argument_parser_test", "-i", "abc"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // fail on number followed by character
    char const * argv2[] = {"./argument_parser_test", "-i", "2abc"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);

    // fail on double
    char const * argv3[] = {"./argument_parser_test", "-i", "3.12"};
    seqan3::argument_parser parser3{"test_parser", 3, argv3, seqan3::update_notifications::off};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser3.parse(), seqan3::argument_parser_error);

    // fail on negative number for unsigned
    unsigned option_value_u;
    char const * argv4[] = {"./argument_parser_test", "-i", "-1"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, seqan3::update_notifications::off};
    parser4.add_option(option_value_u, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser4.parse(), seqan3::argument_parser_error);

    // fail on overflow
    int8_t option_value_int8;
    char const * argv5[] = {"./argument_parser_test", "-i", "129"};
    seqan3::argument_parser parser5{"test_parser", 3, argv5, seqan3::update_notifications::off};
    parser5.add_option(option_value_int8, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser5.parse(), seqan3::argument_parser_error);

    uint8_t option_value_uint8;
    char const * argv6[] = {"./argument_parser_test", "-i", "267"};
    seqan3::argument_parser parser6{"test_parser", 3, argv6, seqan3::update_notifications::off};
    parser6.add_option(option_value_uint8, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser6.parse(), seqan3::argument_parser_error);
}

TEST(parse_type_test, parse_error_double_option)
{
    double option_value;

    // fail on character
    char const * argv[] = {"./argument_parser_test", "-i", "abc"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'd', "double-option", "this is a double option.");

    EXPECT_THROW(parser.parse(), seqan3::argument_parser_error);

    // fail on number followed by character
    char const * argv2[] = {"./argument_parser_test", "-d", "12.457a"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'd', "double-option", "this is a double option.");

    EXPECT_THROW(parser2.parse(), seqan3::argument_parser_error);
}

TEST(parse_test, too_many_arguments_error)
{
    int option_value;

    char const * argv[] = {"./argument_parser_test", "5", "15"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_positional_option(option_value, "this is an int option.");

    EXPECT_THROW(parser.parse(), seqan3::too_many_arguments);

    // since -- indicates -i as a positional argument, this causes a too many args error
    char const * argv2[] = {"./argument_parser_test", "2", "--", "-i"};
    seqan3::argument_parser parser2{"test_parser", 4, argv2, seqan3::update_notifications::off};
    parser2.add_positional_option(option_value, "normal int positional argument.");
    parser2.add_option(option_value, 'i', "int-option", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::too_many_arguments);
}

TEST(parse_test, too_few_arguments_error)
{
    int option_value;

    char const * argv[] = {"./argument_parser_test", "15"};
    seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
    parser.add_positional_option(option_value, "this is an int option.");
    parser.add_positional_option(option_value, "this is another option.");

    EXPECT_THROW(parser.parse(), seqan3::too_few_arguments);

    // since -- indicates -i as a positional argument, this causes a too many args error
    char const * argv2[] = {"./argument_parser_test", "-i", "2"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    parser2.add_positional_option(option_value, "normal int positional argument.");
    parser2.add_option(option_value, 'i', "int-option", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::too_few_arguments);
}

TEST(parse_test, unknown_option_error)
{
    // unknown short option
    char const * argv[] = {"./argument_parser_test", "-i", "15"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};

    EXPECT_THROW(parser.parse(), seqan3::unknown_option);

    // unknown long option
    char const * argv2[] = {"./argument_parser_test", "--arg", "8"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};

    EXPECT_THROW(parser2.parse(), seqan3::unknown_option);

    // unknown short flag
    char const * argv3[] = {"./argument_parser_test", "-a"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};

    EXPECT_THROW(parser3.parse(), seqan3::unknown_option);

    // unknown long flag
    char const * argv4[] = {"./argument_parser_test", "--arg"};
    seqan3::argument_parser parser4{"test_parser", 2, argv4, seqan3::update_notifications::off};

    EXPECT_THROW(parser4.parse(), seqan3::unknown_option);

    // negative numbers are seen as options
    char const * argv5[] = {"./argument_parser_test", "-5"};
    seqan3::argument_parser parser5{"test_parser", 2, argv5, seqan3::update_notifications::off};

    EXPECT_THROW(parser5.parse(), seqan3::unknown_option);

    // unknown short option in more complex command line
    int option_value_i;
    std::string option_value_a;
    std::string positional_option;
    char const * argv6[] = {"./argument_parser_test", "-i", "129", "arg1", "-b", "bcd", "-a", "abc"};
    seqan3::argument_parser parser6{"test_parser", 8, argv6, seqan3::update_notifications::off};
    parser6.add_option(option_value_i, 'i', "int-option", "this is a int option.");
    parser6.add_option(option_value_a, 'a', "string-option", "this is a string option.");
    parser6.add_positional_option(positional_option, "normal int positional argument.");

    EXPECT_THROW(parser6.parse(), seqan3::unknown_option);
}

TEST(parse_test, option_declared_multiple_times_error)
{
    int option_value;

    // short option
    char const * argv[] = {"./argument_parser_test", "-i", "15", "-i", "3"};
    seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'i', "long", "this is a int option.");

    EXPECT_THROW(parser.parse(), seqan3::option_declared_multiple_times);

    // since -- indicates -i as a positional argument, this causes a too many args error
    char const * argv2[] = {"./argument_parser_test", "--long", "5", "--long", "6"};
    seqan3::argument_parser parser2{"test_parser", 5, argv2, seqan3::update_notifications::off};
    parser2.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser2.parse(), seqan3::option_declared_multiple_times);

    // since -- indicates -i as a positional argument, this causes a too many args error
    char const * argv3[] = {"./argument_parser_test", "-i", "5", "--long", "6"};
    seqan3::argument_parser parser3{"test_parser", 5, argv3, seqan3::update_notifications::off};
    parser3.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser3.parse(), seqan3::option_declared_multiple_times);
}

TEST(parse_test, required_option_missing)
{
    int option_value;

    // option is required
    char const * argv[] = {"./argument_parser_test", "5", "-i", "15"};
    seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'i', "int-option", "this is an int option.");
    parser.add_option(option_value, 'a', "req-option", "I am required.", seqan3::option_spec::required);
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
    seqan3::argument_parser parser{"test_parser", 2, argv_all_const, seqan3::update_notifications::off};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // none const
    flag_value = false;
    parser = seqan3::argument_parser{"test_parser", 2, argv, seqan3::update_notifications::off};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 1
    flag_value = false;
    char const * argv_const1[] = {"./argument_parser_test", "-f"};
    parser = seqan3::argument_parser{"test_parser", 2, argv_const1, seqan3::update_notifications::off};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 2
    flag_value = false;
    char * const argv_const2[] = {arg1, arg2};
    parser = seqan3::argument_parser{"test_parser", 2, argv_const2, seqan3::update_notifications::off};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 12
    flag_value = false;
    char const * const argv_const12[] = {arg1, arg2};
    parser = seqan3::argument_parser{"test_parser", 2, argv_const12, seqan3::update_notifications::off};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);
}

TEST(parse_test, multiple_empty_options)
{
    int option_value;

    {
        char const * argv[]{"./empty_long", "-s=1"};
        seqan3::argument_parser parser{"empty_long", 2, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'i', "", "no long");
        parser.add_option(option_value, 's', "", "no long");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(1, option_value);
    }

    {
        char const * argv[]{"./empty_long", "-s=1", "--unknown"};
        seqan3::argument_parser parser{"empty_long", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'i', "", "no long");
        parser.add_option(option_value, 's', "", "no long");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    {
        char const * argv[]{"./empty_short", "--long=2"};
        seqan3::argument_parser parser{"empty_short", 2, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, '\0', "longi", "no short");
        parser.add_option(option_value, '\0', "long", "no short");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(2, option_value);
    }
}

TEST(parse_test, version_check_option_error)
{
    { // version-check must be followed by a value
        char const * argv[] = {"./argument_parser_test", "--version-check"};
        EXPECT_THROW((seqan3::argument_parser{"test_parser", 2, argv}), seqan3::argument_parser_error);
    }

    { // version-check value must be 0 or 1
        char const * argv[] = {"./argument_parser_test", "--version-check", "foo"};
        EXPECT_THROW((seqan3::argument_parser{"test_parser", 3, argv}), seqan3::argument_parser_error);
    }
}

TEST(parse_test, subcommand_argument_parser_success)
{
    bool flag_value{false};
    std::string option_value{};

    // parsing
    {
        char const * argv[]{"./top_level", "-f", "sub1", "foo"};
        seqan3::argument_parser top_level_parser{"top_level",
                                                 4,
                                                 argv,
                                                 seqan3::update_notifications::off,
                                                 {"sub1", "sub2"}};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        EXPECT_NO_THROW(top_level_parser.parse());
        EXPECT_EQ(true, flag_value);

        seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser();

        EXPECT_EQ(sub_parser.info.app_name, "top_level-sub1");

        sub_parser.add_positional_option(option_value, "foo bar");

        EXPECT_NO_THROW(sub_parser.parse());
        EXPECT_EQ("foo", option_value);
    }

    flag_value = false; // reinstate to default value

    // top-level help page
    {
        char const * argv[]{"./top_level", "-h", "-f", "sub1", "foo"};
        seqan3::argument_parser top_level_parser{"top_level",
                                                 5,
                                                 argv,
                                                 seqan3::update_notifications::off,
                                                 {"sub1", "sub2"}};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        testing::internal::CaptureStdout();
        EXPECT_EXIT(top_level_parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        EXPECT_FALSE(std::string{testing::internal::GetCapturedStdout()}.empty());
    }

    flag_value = false; // reinstate to default value

    // sub-parser help page
    {
        char const * argv[]{"./top_level", "-f", "sub1", "-h"};
        seqan3::argument_parser top_level_parser{"top_level",
                                                 4,
                                                 argv,
                                                 seqan3::update_notifications::off,
                                                 {"sub1", "sub2"}};
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
    char const * argv[]{"./top_level", "subiddysub", "-f"};
    { // see issue https://github.com/seqan/seqan3/issues/2172
        seqan3::argument_parser top_level_parser{"top_level",
                                                 3,
                                                 argv,
                                                 seqan3::update_notifications::off,
                                                 {"sub1", "sub2"}};

        EXPECT_THROW(top_level_parser.parse(), seqan3::argument_parser_error);
    }

    // sub command can contain dash, see https://github.com/seqan/product_backlog/issues/234
    {
        EXPECT_NO_THROW((seqan3::argument_parser{"top_level", 2, argv, seqan3::update_notifications::off, {"-dash"}}));
    }
}

TEST(parse_test, issue1544)
{
    { // wrong separation of long value:
        std::string option_value;
        char const * argv[] = {"./argument_parser_test", "--foohallo"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'f', "foo", "this is a string option.");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    { // unknown option (`--foo-bar`) that has a prefix of a known option (`--foo`)
        std::string option_value;
        char const * argv[] = {"./argument_parser_test", "--foo", "hallo", "--foo-bar", "ballo"};
        seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'f', "foo", "this is a string option.");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    { // known option (`--foo-bar`) that has a prefix of a unknown option (`--foo`)
        std::string option_value;
        char const * argv[] = {"./argument_parser_test", "--foo", "hallo", "--foo-bar", "ballo"};
        seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'f', "foo-bar", "this is a string option.");

        EXPECT_THROW(parser.parse(), seqan3::unknown_option);
    }

    { // known option (`--foo`) is a prefix of another known option (`--foo-bar`)
        std::string foo_option_value;
        std::string foobar_option_value;
        char const * argv[] = {"./argument_parser_test", "--foo", "hallo", "--foo-bar", "ballo"};
        seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
        parser.add_option(foo_option_value, 'f', "foo", "this is a prefix of foobar.");
        parser.add_option(foobar_option_value, 'b', "foo-bar", "this has prefix foo.");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(foo_option_value, "hallo");
        EXPECT_EQ(foobar_option_value, "ballo");
    }
}

TEST(parse_test, is_option_set)
{
    std::string option_value{};
    char const * argv[] = {"./argument_parser_test", "-l", "hallo", "--foobar", "ballo", "--", "--loo"};
    seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'l', "loo", "this is a option.");
    parser.add_option(option_value, 'f', "foobar", "this is a option.");

    EXPECT_THROW(parser.is_option_set("foo"), seqan3::design_error); // you cannot call option_is_set before parse()

    EXPECT_NO_THROW(parser.parse());

    EXPECT_TRUE(parser.is_option_set('l'));
    EXPECT_TRUE(parser.is_option_set("foobar"));

    EXPECT_FALSE(parser.is_option_set('f'));
    EXPECT_FALSE(parser.is_option_set("loo")); // --loo is behind the `--` which signals the end of options!

    // errors:
    EXPECT_THROW(parser.is_option_set("l"), seqan3::design_error); // short identifiers are passed as chars not strings
    EXPECT_THROW(parser.is_option_set("f"), seqan3::design_error); // short identifiers are passed as chars not strings

    EXPECT_THROW(parser.is_option_set("foo"), seqan3::design_error);
    EXPECT_THROW(parser.is_option_set("--"), seqan3::design_error);
    EXPECT_THROW(parser.is_option_set(""), seqan3::design_error);

    EXPECT_THROW(parser.is_option_set('!'), seqan3::design_error);
    EXPECT_THROW(parser.is_option_set('-'), seqan3::design_error);
    EXPECT_THROW(parser.is_option_set('_'), seqan3::design_error);
    EXPECT_THROW(parser.is_option_set('\0'), seqan3::design_error);
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
    static inline std::unordered_map<std::string_view, Other::bar> const enumeration_names{{"one", Other::bar::one},
                                                                                           {"1", Other::bar::one},
                                                                                           {"two", Other::bar::two},
                                                                                           {"2", Other::bar::two}};
};
} // namespace seqan3::custom

TEST(parse_type_test, parse_success_enum_option)
{
    {
        foo::bar option_value{};

        char const * argv[] = {"./argument_parser_test", "-e", "two"};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE(option_value == foo::bar::two);
    }

    {
        Other::bar option_value{};

        char const * argv[] = {"./argument_parser_test", "-e", "two"};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE(option_value == Other::bar::two);
    }
}

TEST(parse_type_test, parse_error_enum_option)
{
    foo::bar option_value{};

    char const * argv[] = {"./argument_parser_test", "-e", "four"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

    EXPECT_THROW(parser.parse(), seqan3::user_input_error);
}

// https://github.com/seqan/seqan3/issues/2464
TEST(parse_test, issue2464)
{
    using option_t = foo::bar;
    // Using a non-existing value of foo::bar should throw.
    {
        char const * argv[] = {"./argument_parser_test", "-e", "nine"};

        option_t option_value{};

        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");
        EXPECT_THROW(parser.parse(), seqan3::user_input_error);
    }
    {
        char const * argv[] = {"./argument_parser_test", "-e", "one", "-e", "nine"};

        std::vector<option_t> option_values{};

        seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'e', "enum-option", "this is an enum option.");
        EXPECT_THROW(parser.parse(), seqan3::user_input_error);
    }

    // Invalid inputs for enums are handled before any validator is evaluated.
    // Thus the exception will be seqan3::user_input_error and not seqan3::validation_error.
    {
        char const * argv[] = {"./argument_parser_test", "-e", "nine"};

        seqan3::value_list_validator enum_validator{(seqan3::enumeration_names<option_t> | std::views::values)};
        option_t option_value{};

        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value,
                          'e',
                          "enum-option",
                          "this is an enum option.",
                          seqan3::option_spec::advanced,
                          enum_validator);
        EXPECT_THROW(parser.parse(), seqan3::user_input_error);
    }
    {
        char const * argv[] = {"./argument_parser_test", "-e", "one", "-e", "nine"};

        seqan3::value_list_validator enum_validator{(seqan3::enumeration_names<option_t> | std::views::values)};
        std::vector<option_t> option_values{};

        seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
        parser.add_option(option_values,
                          'e',
                          "enum-option",
                          "this is an enum option.",
                          seqan3::option_spec::advanced,
                          enum_validator);
        EXPECT_THROW(parser.parse(), seqan3::user_input_error);
    }
}

TEST(parse_test, enum_error_message)
{
    // foo::bar does not contain duplicate values
    {
        char const * argv[] = {"./argument_parser_test", "-e", "nine"};

        foo::bar option_value{};

        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

        std::string expected_message{"You have chosen an invalid input value: nine. "
                                     "Please use one of: [one,two,three]"};

        try
        {
            parser.parse();
            FAIL();
        }
        catch (seqan3::user_input_error const & exception)
        {
            EXPECT_EQ(expected_message, exception.what());
        }
        catch (...)
        {
            FAIL();
        }
    }
    // Other::bar does contain duplicate values
    {
        char const * argv[] = {"./argument_parser_test", "-e", "nine"};

        Other::bar option_value{};

        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_value, 'e', "enum-option", "this is an enum option.");

        std::string expected_message{"You have chosen an invalid input value: nine. "
                                     "Please use one of: [1,one,2,two]"};

        try
        {
            parser.parse();
            FAIL();
        }
        catch (seqan3::user_input_error const & exception)
        {
            EXPECT_EQ(expected_message, exception.what());
        }
        catch (...)
        {
            FAIL();
        }
    }
}

// https://github.com/seqan/seqan3/issues/2835
TEST(parse_test, error_message_parsing)
{
    char const * argv[] = {"./argument_parser_test", "--value", "-30"};

    uint64_t option_value{};

    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    parser.add_option(option_value, '\0', "value", "Please specify a value.");

    std::string expected_message{"Value parse failed for --value: Argument -30 could not be parsed as type "
                                 "unsigned 64 bit integer."};

    try
    {
        parser.parse();
        FAIL();
    }
    catch (seqan3::user_input_error const & exception)
    {
        EXPECT_EQ(expected_message, exception.what());
    }
    catch (...)
    {
        FAIL();
    }
}

// https://github.com/seqan/seqan3/pull/2381
TEST(parse_test, container_options)
{
    {
        std::vector<foo::bar> option_values{};

        char const * argv[] = {"./argument_parser_test", "-e", "two", "-e", "one", "-e", "three"};
        seqan3::argument_parser parser{"test_parser", 7, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'e', "enum-option", "this is an enum option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<foo::bar>{foo::bar::two, foo::bar::one, foo::bar::three}));
    }

    {
        std::vector<int> option_values{};

        char const * argv[] = {"./argument_parser_test", "-i", "2", "-i", "1", "-i", "3"};
        seqan3::argument_parser parser{"test_parser", 7, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'i', "int-option", "this is an int option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<int>{2, 1, 3}));
    }

    {
        std::vector<bool> option_values{};

        char const * argv[] = {"./argument_parser_test", "-b", "true", "-b", "false", "-b", "true"};
        seqan3::argument_parser parser{"test_parser", 7, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'b', "bool-option", "this is a bool option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<bool>{true, false, true}));
    }
}

// https://github.com/seqan/seqan3/issues/2393
TEST(parse_test, container_default)
{
    // overwrite default
    {
        std::vector<int> option_values{1, 2, 3};

        char const * argv[] = {"./argument_parser_test", "-i", "2", "-i", "1", "-i", "3"};
        seqan3::argument_parser parser{"test_parser", 7, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'i', "int-option", "this is an int option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<int>{2, 1, 3}));
    }
    // overwrite default, parameters are not consecutive
    {
        std::vector<int> option_values{1, 2, 3};
        bool bool_opt{false};

        char const * argv[] = {"./argument_parser_test", "-i", "2", "-b", "true", "-i", "1", "-i", "3"};
        seqan3::argument_parser parser{"test_parser", 9, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'i', "int-option", "this is an int option.");
        parser.add_option(bool_opt, 'b', "bool-option", "this is a bool option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<int>{2, 1, 3}));
    }
    // use default
    {
        std::vector<int> option_values{1, 2, 3};
        bool bool_opt{false};

        char const * argv[] = {"./argument_parser_test", "-b", "true"};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        parser.add_option(option_values, 'i', "int-option", "this is an int option.");
        parser.add_option(bool_opt, 'b', "bool-option", "this is a bool option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<int>{1, 2, 3}));
    }
    // overwrite default for positional options
    {
        std::vector<int> option_values{1, 2, 3};

        char const * argv[] = {"./argument_parser_test", "2", "1", "3"};
        seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
        parser.add_positional_option(option_values, "this is an int option.");

        EXPECT_NO_THROW(parser.parse());

        EXPECT_TRUE(option_values == (std::vector<int>{2, 1, 3}));
    }
}
