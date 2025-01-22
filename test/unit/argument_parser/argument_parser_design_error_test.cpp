// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/argument_parser/argument_parser.hpp>

TEST(design_error, app_name_validation)
{
    char const * argv[] = {"./argument_parser_test"};

    EXPECT_NO_THROW((seqan3::argument_parser{"test_parser", 1, argv}));
    EXPECT_NO_THROW((seqan3::argument_parser{"test-parser1234_foo", 1, argv}));

    EXPECT_THROW((seqan3::argument_parser{"test parser", 1, argv}), seqan3::design_error);
    EXPECT_THROW((seqan3::argument_parser{"test;", 1, argv}), seqan3::design_error);
    EXPECT_THROW((seqan3::argument_parser{";", 1, argv}), seqan3::design_error);
    EXPECT_THROW((seqan3::argument_parser{"test;bad script:D", 1, argv}), seqan3::design_error);
}

TEST(parse_test, design_error)
{
    int option_value;

    // short option
    char const * argv[] = {"./argument_parser_test"};
    seqan3::argument_parser parser{"test_parser", 1, argv};
    parser.add_option(option_value, 'i', "int", "this is a int option.");
    EXPECT_THROW(parser.add_option(option_value, 'i', "aint", "oh oh same id."), seqan3::design_error);

    // long option
    seqan3::argument_parser parser2{"test_parser", 1, argv};
    parser2.add_option(option_value, 'i', "int", "this is an int option.");
    EXPECT_THROW(parser2.add_option(option_value, 'a', "int", "oh oh another id."), seqan3::design_error);

    // empty identifier
    seqan3::argument_parser parser3{"test_parser", 1, argv};
    EXPECT_THROW(parser3.add_option(option_value, '\0', "", "oh oh all is empty."), seqan3::design_error);

    bool true_value{true};

    // default true
    seqan3::argument_parser parser4{"test_parser", 1, argv};
    EXPECT_THROW(parser4.add_flag(true_value, 'i', "int", "oh oh default is true."), seqan3::design_error);

    bool flag_value{false};

    // short flag
    seqan3::argument_parser parser5{"test_parser", 1, argv};
    parser5.add_flag(flag_value, 'i', "int1", "this is an int option.");
    EXPECT_THROW(parser5.add_flag(flag_value, 'i', "int2", "oh oh another id."), seqan3::design_error);

    // long flag
    seqan3::argument_parser parser6{"test_parser", 1, argv};
    parser6.add_flag(flag_value, 'i', "int", "this is an int option.");
    EXPECT_THROW(parser6.add_flag(flag_value, 'a', "int", "oh oh another id."), seqan3::design_error);

    // empty identifier
    seqan3::argument_parser parser7{"test_parser", 1, argv};
    EXPECT_THROW(parser7.add_flag(flag_value, '\0', "", "oh oh another id."), seqan3::design_error);

    // positional option not at the end
    char const * argv2[] = {"./argument_parser_test", "arg1", "arg2", "arg3"};
    std::vector<int> vec;
    seqan3::argument_parser parser8{"test_parser", 4, argv2};
    parser8.add_positional_option(vec, "oh oh list not at the end.");
    EXPECT_THROW(parser8.add_positional_option(option_value, "desc."), seqan3::design_error);

    // using h, help, advanced-help, and export-help
    seqan3::argument_parser parser9{"test_parser", 1, argv};
    EXPECT_THROW(parser9.add_option(option_value, 'h', "", "-h is bad."), seqan3::design_error);
    EXPECT_THROW(parser9.add_option(option_value, '\0', "help", "help is bad."), seqan3::design_error);
    EXPECT_THROW(parser9.add_option(option_value, '\0', "advanced-help", "advanced-help is bad"), seqan3::design_error);
    EXPECT_THROW(parser9.add_option(option_value, '\0', "export-help", "export-help is bad"), seqan3::design_error);

    // using one-letter long identifiers.
    seqan3::argument_parser parser10{"test_parser", 1, argv};
    EXPECT_THROW(parser10.add_option(option_value, 'y', "z", "long identifier is one letter"), seqan3::design_error);
    EXPECT_THROW(parser10.add_flag(flag_value, 'y', "z", "long identifier is one letter"), seqan3::design_error);

    // using non-printable characters
    seqan3::argument_parser parser11{"test_parser", 1, argv};
    EXPECT_THROW(parser11.add_option(option_value, '\t', "no\n", "tab and newline don't work!"), seqan3::design_error);
    EXPECT_THROW(parser11.add_flag(flag_value, 'i', "no\n", "tab and newline don't work!"), seqan3::design_error);
    EXPECT_THROW(parser11.add_flag(flag_value, 'a', "-no", "can't start long_id with a hyphen"), seqan3::design_error);
}

TEST(parse_test, parse_called_twice)
{
    std::string option_value;

    char const * argv[] = {"./argument_parser_test", "--version-check", "false", "-s", "option_string"};
    seqan3::argument_parser parser{"test_parser", 5, argv};
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    EXPECT_THROW(parser.parse(), seqan3::design_error);
}

TEST(parse_test, subcommand_argument_parser_error)
{
    bool flag_value{};

    // subcommand parsing was not enabled on construction but get_sub_parser() is called
    {
        char const * argv[]{"./top_level", "-f"};
        seqan3::argument_parser top_level_parser{"top_level", 2, argv, seqan3::update_notifications::off};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        EXPECT_NO_THROW(top_level_parser.parse());
        EXPECT_EQ(true, flag_value);

        EXPECT_THROW(top_level_parser.get_sub_parser(), seqan3::design_error);
    }

    // subcommand key word must only contain alpha numeric characters
    {
        char const * argv[]{"./top_level", "-f"};
        EXPECT_THROW((seqan3::argument_parser{"top_level", 2, argv, seqan3::update_notifications::off, {"with space"}}),
                     seqan3::design_error);
    }

    // no positional/options are allowed
    {
        char const * argv[]{"./top_level", "foo"};
        seqan3::argument_parser top_level_parser{"top_level", 2, argv, seqan3::update_notifications::off, {"foo"}};

        EXPECT_THROW((top_level_parser.add_option(flag_value, 'f', "foo", "foo bar")), seqan3::design_error);
        EXPECT_THROW((top_level_parser.add_positional_option(flag_value, "foo bar")), seqan3::design_error);
    }
}
