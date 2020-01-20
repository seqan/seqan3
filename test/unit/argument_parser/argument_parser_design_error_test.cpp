// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

TEST(parser_design_error, app_name_validation)
{
    const char * argv[] = {"./argument_parser_test"};

    EXPECT_NO_THROW((argument_parser{"test_parser", 1, argv}));
    EXPECT_NO_THROW((argument_parser{"test-parser1234_foo", 1, argv}));

    EXPECT_THROW((argument_parser{"test parser", 1, argv}),       parser_design_error);
    EXPECT_THROW((argument_parser{"test;", 1, argv}),             parser_design_error);
    EXPECT_THROW((argument_parser{";", 1, argv}),                 parser_design_error);
    EXPECT_THROW((argument_parser{"test;bad script:D", 1, argv}), parser_design_error);
}

TEST(parse_test, parser_design_error)
{
    int option_value;

    // short option
    const char * argv[] = {"./argument_parser_test"};
    argument_parser parser{"test_parser", 1, argv};
    parser.add_option(option_value, 'i', "int", "this is a int option.");
    EXPECT_THROW(parser.add_option(option_value, 'i', "aint", "oh oh same id."),
                 parser_design_error);

    // long option
    argument_parser parser2{"test_parser", 1, argv};
    parser2.add_option(option_value, 'i', "int", "this is an int option.");
    EXPECT_THROW(parser2.add_option(option_value, 'a', "int", "oh oh another id."),
                 parser_design_error);

    // empty identifier
    argument_parser parser3{"test_parser", 1, argv};
    EXPECT_THROW(parser3.add_option(option_value, '\0', "", "oh oh all is empty."),
                 parser_design_error);

    bool flag_value;

    // short flag
    argument_parser parser4{"test_parser", 1, argv};
    parser4.add_flag(flag_value, 'i', "int1", "this is an int option.");
    EXPECT_THROW(parser4.add_flag(flag_value, 'i', "int2", "oh oh another id."),
                 parser_design_error);

    // long flag
    argument_parser parser5{"test_parser", 1, argv};
    parser5.add_flag(flag_value, 'i', "int", "this is an int option.");
    EXPECT_THROW(parser5.add_flag(flag_value, 'a', "int", "oh oh another id."),
                 parser_design_error);

    // empty identifier
    argument_parser parser6{"test_parser", 1, argv};
    EXPECT_THROW(parser6.add_flag(flag_value, '\0', "", "oh oh another id."),
                 parser_design_error);

    // positional option not at the end
    const char * argv2[] = {"./argument_parser_test", "arg1", "arg2", "arg3"};
    std::vector<int> vec;
    argument_parser parser7{"test_parser", 4, argv2};
    parser7.add_positional_option(vec, "oh oh list not at the end.");
    EXPECT_THROW(parser7.add_positional_option(option_value, "desc."), parser_design_error);

    // using h, help, advanced-help, and export-help
    argument_parser parser8{"test_parser", 1, argv};
    EXPECT_THROW(parser8.add_option(option_value, 'h', "", "-h is bad."),
                 parser_design_error);
    EXPECT_THROW(parser8.add_option(option_value, '\0', "help", "help is bad."),
                 parser_design_error);
    EXPECT_THROW(parser8.add_option(option_value, '\0', "advanced-help",
                 "advanced-help is bad"), parser_design_error);
    EXPECT_THROW(parser8.add_option(option_value, '\0', "export-help",
                 "export-help is bad"), parser_design_error);

    // using one-letter long identifiers.
    argument_parser parser9{"test_parser", 1, argv};
    EXPECT_THROW(parser9.add_option(option_value, 'y', "z", "long identifier is one letter"),
                 parser_design_error);
    EXPECT_THROW(parser9.add_flag(flag_value, 'y', "z", "long identifier is one letter"),
                 parser_design_error);

    // using non-printable characters
    argument_parser parser10{"test_parser", 1, argv};
    EXPECT_THROW(parser10.add_option(option_value, '\t', "no\n", "tab and newline don't work!"),
                 parser_design_error);
    EXPECT_THROW(parser10.add_flag(flag_value, 'i', "no\n", "tab and newline don't work!"),
                 parser_design_error);
    EXPECT_THROW(parser10.add_flag(flag_value, 'a', "-no", "can't start long_id with a hyphen"),
                 parser_design_error);
}

TEST(parse_test, parse_called_twice)
{
    std::string option_value;

    const char * argv[] = {"./argument_parser_test", "--version-check", "0", "-s", "option_string"};
    argument_parser parser{"test_parser", 5, argv};
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    EXPECT_THROW(parser.parse(), parser_design_error);
}

TEST(parse_test, subcommand_argument_parser_error)
{
    bool flag_value{};

    // subcommand parsing was not enabled on construction but get_sub_parser() is called
    {
        const char * argv[]{"./top_level", "-f"};
        argument_parser top_level_parser{"top_level", 2, argv, false};
        top_level_parser.add_flag(flag_value, 'f', "foo", "foo bar");

        EXPECT_NO_THROW(top_level_parser.parse());
        EXPECT_EQ(true, flag_value);

        EXPECT_THROW(top_level_parser.get_sub_parser(), parser_design_error);
    }

    // subcommand key word must only contain alpha numeric characters
    {
        const char * argv[]{"./top_level", "-f"};
        EXPECT_THROW((argument_parser{"top_level", 2, argv, false, {"with space"}}), parser_design_error);
        EXPECT_THROW((argument_parser{"top_level", 2, argv, false, {"-dash"}}), parser_design_error);
    }

    // no positional/options are allowed
    {
        const char * argv[]{"./top_level", "foo"};
        argument_parser top_level_parser{"top_level", 2, argv, false, {"foo"}};

        EXPECT_THROW((top_level_parser.add_option(flag_value, 'f', "foo", "foo bar")), parser_design_error);
        EXPECT_THROW((top_level_parser.add_positional_option(flag_value, "foo bar")), parser_design_error);
    }
}
