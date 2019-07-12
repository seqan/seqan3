// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

TEST(parse_type_test, add_option_short_id)
{
    std::string option_value;

    const char * argv[] = {"./argument_parser_test", "-s", "option_string"};
    argument_parser parser{"test_parser", 3, argv, false};
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
    argument_parser parser2{"test_parser", 2, argv2, false};
    parser2.add_option(option_value, 'S', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    // ad with = sign
    const char * argv3[] = {"./argument_parser_test", "-s=option_string"};
    argument_parser parser3{"test_parser", 2, argv3, false};
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
    argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    // add with no space
    const char * argv2[] = {"./argument_parser_test", "--string-optionoption_string"};
    argument_parser parser2{"test_parser", 2, argv2, false};
    parser2.add_option(option_value, 'S', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    const char * argv3[] = {"./argument_parser_test", "--string-option=option_string"};
    argument_parser parser3{"test_parser", 2, argv3, false};
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
    argument_parser parser{"test_parser", 2, argv, false};
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
    argument_parser parser{"test_parser", 2, argv, false};
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
    argument_parser parser{"test_parser", 2, argv, false};
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
    argument_parser parser{"test_parser", 2, argv, false};
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
    argument_parser parser{"test_parser", 5, argv, false};
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
    argument_parser parser2{"test_parser", 5, argv, false};
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
    argument_parser parser3{"test_parser", 5, argv, false};
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
    argument_parser parser4{"test_parser", 5, argv, false};
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
    argument_parser parser5{"test_parser", 5, argv, false};
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
    argument_parser parser6{"test_parser", 5, argv, false};
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
    argument_parser parser{"test_parser", 5, argv, false};
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
    argument_parser parser2{"test_parser", 5, argv2, false};
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
    argument_parser parser3{"test_parser", 5, argv3, false};
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
    argument_parser parser4{"test_parser", 5, argv4, false};
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
    argument_parser parser5{"test_parser", 5, argv5, false};
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
    argument_parser parser6{"test_parser", 5, argv6, false};
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
    argument_parser parser{"test_parser", 3, argv, false};
    parser.add_positional_option(option_value, "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "-strange");

    // negative integer option
    int option_value_int;
    const char * argv2[] = {"./argument_parser_test", "--", "-120"};
    argument_parser parser2{"test_parser", 3, argv2, false};
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
    argument_parser parser{"test_parser", 3, argv, false};
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
    argument_parser parser{"test_parser", 2, argv, false};
    parser.add_option(option_value, 'i', "long", "this is a int option.");

    EXPECT_THROW(parser.parse(), parser_invalid_argument);

    // long option
    const char * argv2[] = {"./argument_parser_test", "--long"};
    argument_parser parser2{"test_parser", 2, argv2, false};
    parser2.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser2.parse(), parser_invalid_argument);

    // short option
    const char * argv3[] = {"./argument_parser_test", "-i="};
    argument_parser parser3{"test_parser", 2, argv3, false};
    parser3.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser3.parse(), parser_invalid_argument);

    // short option
    const char * argv4[] = {"./argument_parser_test", "--long="};
    argument_parser parser4{"test_parser", 2, argv4, false};
    parser4.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser4.parse(), parser_invalid_argument);
}

TEST(parse_type_test, parse_success_bool_option)
{
    bool option_value;
    bool positional_value;

    // numbers 0 and 1
    {
        const char * argv[] = {"./argument_parser_test", "-b", "1", "0"};
        argument_parser parser{"test_parser", 4, argv, false};
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
        argument_parser parser{"test_parser", 4, argv, false};
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
    argument_parser parser{"test_parser", 4, argv, false};
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
    argument_parser parser{"test_parser", 4, argv, false};
    parser.add_option(option_value, 'd', "double-option", "this is a double option.");
    parser.add_positional_option(positional_value, "this is a double positional.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());

    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_FLOAT_EQ(option_value, 12.457);
    EXPECT_FLOAT_EQ(positional_value, 0.123);

    // double expression with e
    const char * argv2[] = {"./argument_parser_test", "-d", "6.0221418e23"};
    argument_parser parser2{"test_parser", 3, argv2, false};
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
    argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'b', "bool-option", "this is a bool option.");

    EXPECT_THROW(parser.parse(), parser_invalid_argument);

    // fail on number input expect 0 and 1
    const char * argv2[] = {"./argument_parser_test", "-b", "124"};
    argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'b', "bool-option", "this is a bool option.");

    EXPECT_THROW(parser2.parse(), parser_invalid_argument);
}

TEST(parse_type_test, parse_error_int_option)
{
    int option_value;

    // fail on character
    const char * argv[] = {"./argument_parser_test", "-i", "abc"};
    argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser.parse(), parser_invalid_argument);

    // fail on number followed by character
    const char * argv2[] = {"./argument_parser_test", "-i", "2abc"};
    argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser2.parse(), parser_invalid_argument);

    // fail on double
    const char * argv3[] = {"./argument_parser_test", "-i", "3.12"};
    argument_parser parser3{"test_parser", 3, argv3, false};
    parser3.add_option(option_value, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser3.parse(), parser_invalid_argument);

    // fail on negative number for unsigned
    unsigned option_value_u;
    const char * argv4[] = {"./argument_parser_test", "-i", "-1"};
    argument_parser parser4{"test_parser", 3, argv4, false};
    parser4.add_option(option_value_u, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser4.parse(), parser_invalid_argument);

    // fail on overflow
    int8_t option_value_int8;
    const char * argv5[] = {"./argument_parser_test", "-i", "129"};
    argument_parser parser5{"test_parser", 3, argv5, false};
    parser5.add_option(option_value_int8, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser5.parse(), parser_invalid_argument);

    uint8_t option_value_uint8;
    const char * argv6[] = {"./argument_parser_test", "-i", "267"};
    argument_parser parser6{"test_parser", 3, argv6, false};
    parser6.add_option(option_value_uint8, 'i', "int-option", "this is a int option.");

    EXPECT_THROW(parser6.parse(), parser_invalid_argument);
}

TEST(parse_type_test, parse_error_double_option)
{
    double option_value;

    // fail on character
    const char * argv[] = {"./argument_parser_test", "-i", "abc"};
    argument_parser parser{"test_parser", 3, argv, false};
    parser.add_option(option_value, 'd', "double-option", "this is a double option.");

    EXPECT_THROW(parser.parse(), parser_invalid_argument);

    // fail on number followed by character
    const char * argv2[] = {"./argument_parser_test", "-d", "12.457a"};
    argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_option(option_value, 'd', "double-option", "this is a double option.");

    EXPECT_THROW(parser2.parse(), parser_invalid_argument);
}


TEST(parse_test, too_many_arguments_error)
{
    int option_value;

    const char * argv[] = {"./argument_parser_test", "5", "15"};
    argument_parser parser{"test_parser", 3, argv, false};
    parser.add_positional_option(option_value, "this is an int option.");

    EXPECT_THROW(parser.parse(), too_many_arguments);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv2[] = {"./argument_parser_test", "2", "--", "-i"};
    argument_parser parser2{"test_parser", 4, argv2, false};
    parser2.add_positional_option(option_value, "normal int positional argument.");
    parser2.add_option(option_value, 'i', "int-option", "this is an int option.");

    EXPECT_THROW(parser2.parse(), too_many_arguments);
}

TEST(parse_test, too_few_arguments_error)
{
    int option_value;

    const char * argv[] = {"./argument_parser_test", "15"};
    argument_parser parser{"test_parser", 2, argv, false};
    parser.add_positional_option(option_value, "this is an int option.");
    parser.add_positional_option(option_value, "this is another option.");

    EXPECT_THROW(parser.parse(), too_few_arguments);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv2[] = {"./argument_parser_test", "-i", "2"};
    argument_parser parser2{"test_parser", 3, argv2, false};
    parser2.add_positional_option(option_value, "normal int positional argument.");
    parser2.add_option(option_value, 'i', "int-option", "this is an int option.");

    EXPECT_THROW(parser2.parse(), too_few_arguments);
}

TEST(parse_test, unknown_option_error)
{
    // unknown short option
    const char * argv[] = {"./argument_parser_test", "-i", "15"};
    argument_parser parser{"test_parser", 3, argv, false};

    EXPECT_THROW(parser.parse(), unknown_option);

    // unknown long option
    const char * argv2[] = {"./argument_parser_test", "--arg", "8"};
    argument_parser parser2{"test_parser", 3, argv2, false};

    EXPECT_THROW(parser2.parse(), unknown_option);

    // unknown short flag
    const char * argv3[] = {"./argument_parser_test", "-a"};
    argument_parser parser3{"test_parser", 2, argv3, false};

    EXPECT_THROW(parser3.parse(), unknown_option);

    // unknown long flag
    const char * argv4[] = {"./argument_parser_test", "--arg"};
    argument_parser parser4{"test_parser", 2, argv4, false};

    EXPECT_THROW(parser4.parse(), unknown_option);

    // negative numbers are seen as options
    const char * argv5[] = {"./argument_parser_test", "-5"};
    argument_parser parser5{"test_parser", 2, argv5, false};

    EXPECT_THROW(parser5.parse(), unknown_option);

    // unknown short option in more complex command line
    int option_value_i;
    std::string option_value_a;
    std::string positional_option;
    const char * argv6[] = {"./argument_parser_test", "-i", "129", "arg1", "-b", "bcd", "-a", "abc"};
    argument_parser parser6{"test_parser", 8, argv6, false};
    parser6.add_option(option_value_i, 'i', "int-option", "this is a int option.");
    parser6.add_option(option_value_a, 'a', "string-option", "this is a string option.");
    parser6.add_positional_option(positional_option, "normal int positional argument.");

    EXPECT_THROW(parser6.parse(), unknown_option);
}

TEST(parse_test, option_declared_multiple_times_error)
{
    int option_value;

    // short option
    const char * argv[] = {"./argument_parser_test", "-i", "15", "-i", "3"};
    argument_parser parser{"test_parser", 5, argv, false};
    parser.add_option(option_value, 'i', "long", "this is a int option.");

    EXPECT_THROW(parser.parse(), option_declared_multiple_times);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv2[] = {"./argument_parser_test", "--long", "5", "--long", "6"};
    argument_parser parser2{"test_parser", 5, argv2, false};
    parser2.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser2.parse(), option_declared_multiple_times);

    // since -- indicates -i as a positional argument, this causes a too many args error
    const char * argv3[] = {"./argument_parser_test", "-i", "5", "--long", "6"};
    argument_parser parser3{"test_parser", 5, argv3, false};
    parser3.add_option(option_value, 'i', "long", "this is an int option.");

    EXPECT_THROW(parser3.parse(), option_declared_multiple_times);
}

TEST(parse_test, required_option_missing)
{
    int option_value;

    // option is required
    const char * argv[] = {"./argument_parser_test", "5", "-i", "15"};
    argument_parser parser{"test_parser", 4, argv, false};
    parser.add_option(option_value, 'i', "int-option", "this is an int option.");
    parser.add_option(option_value, 'a', "req-option", "I am required.", option_spec::REQUIRED);
    parser.add_positional_option(option_value, "this is an int option.");

    EXPECT_THROW(parser.parse(), required_option_missing);
}

TEST(parse_test, argv_const_combinations)
{
    bool flag_value{false};

    char arg1[]{"./argument_parser"};
    char arg2[]{"-f"};
    char * argv[] = {arg1, arg2};

    // all const*
    char const * const * const argv_all_const{argv};
    argument_parser parser{"test_parser", 2, argv_all_const, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // none const
    flag_value = false;
    parser = argument_parser{"test_parser", 2, argv, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 1
    flag_value = false;
    char const * argv_const1[] = {"./argument_parser_test", "-f"};
    parser = argument_parser{"test_parser", 2, argv_const1, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 2
    flag_value = false;
    char * const argv_const2[] = {arg1, arg2};
    parser = argument_parser{"test_parser", 2, argv_const2, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);

    // const 12
    flag_value = false;
    char const * const argv_const12[] = {arg1, arg2};
    parser = argument_parser{"test_parser", 2, argv_const12, false};
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");

    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE(flag_value);
}

TEST(parse_test, multiple_empty_options)
{
    int option_value;

    {
        const char * argv[]{"./empty_long", "-s=1"};
        argument_parser parser{"empty_long", 2, argv, false};
        parser.add_option(option_value, 'i', "", "no long");
        parser.add_option(option_value, 's', "", "no long");

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(1, option_value);
    }

    {
        const char * argv[]{"./empty_long", "-s=1", "--unknown"};
        argument_parser parser{"empty_long", 3, argv, false};
        parser.add_option(option_value, 'i', "", "no long");
        parser.add_option(option_value, 's', "", "no long");

        EXPECT_THROW(parser.parse(), unknown_option);
    }

    {
        const char * argv[]{"./empty_short", "--long=2"};
        argument_parser parser{"empty_short", 2, argv, false};
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
        EXPECT_THROW((argument_parser{"test_parser", 2, argv}), parser_invalid_argument);
    }

    {   // version-check value must be 0 or 1
        const char * argv[] = {"./argument_parser_test", "--version-check", "foo"};
        EXPECT_THROW((argument_parser{"test_parser", 3, argv}), parser_invalid_argument);
    }
}
