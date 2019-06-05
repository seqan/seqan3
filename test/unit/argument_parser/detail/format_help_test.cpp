// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>

#include <gtest/gtest.h>

#include <seqan3/std/ranges>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_help.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/drop.hpp>

using namespace seqan3;

// reused global variables
std::string std_cout;
std::string expected;
int option_value{5};
bool flag_value{};
std::vector<std::string> pos_opt_value{};
const char * argv0[] = {"./help_add_test --version-check 0"};
const char * argv1[] = {"./help_add_test --version-check 0", "-h"};
const char * argv2[] = {"./help_add_test --version-check 0", "-hh"};
const char * argv3[] = {"./help_add_test --version-check 0", "--version"};

std::string const basic_options_str = "OPTIONS"
                                      "Basic options:"
                                      "-h, --help Prints the help page."
                                      "-hh, --advanced-help Prints the help page including advanced options."
                                      "--version Prints the version information."
                                      "--copyright Prints the copyright/license information."
                                      "--export-help (std::string) Export the help page information. "
                                                                   "Value must be one of [html, man]."
                                      "--version-check (bool)"
                                                "Whether to to check for the newest app version. Default: 1.";

std::string const version_str = std::to_string(SEQAN3_VERSION_MAJOR) + "."
                                + std::to_string(SEQAN3_VERSION_MINOR) + "."
                                + std::to_string(SEQAN3_VERSION_PATCH);

std::string const basic_version_str = "VERSION"
                                      "Last update:"
                                      "test_parser version:"
                                      "SeqAn version: " + version_str;

TEST(help_page_printing, short_help)
{
    // Empty call with no options given. For detail::format_short_help
    argument_parser parser0{"empty_options", 1, argv0};
    parser0.info.synopsis.push_back("./some_binary_name synopsis");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser0.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "empty_options"
               "============="
               "./some_binary_name synopsis"
               "Try -h or --help for more information.";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, no_information)
{
    // Empty help call with -h
    argument_parser parser1{"test_parser", 2, argv1};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_options_str +
               basic_version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, with_short_copyright)
{
    // Again, but with short copyright, long copyright, and citation.
    argument_parser short_copy("test_parser", 2, argv1);
    short_copy.info.short_copyright = "short";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(short_copy.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_options_str +
               basic_version_str +
               "LEGAL"
               "test_parser Copyright: short"
               "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, with_long_copyright)
{
    argument_parser long_copy("test_parser", 2, argv1);
    long_copy.info.long_copyright = "long";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(long_copy.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_options_str +
               basic_version_str +
               "LEGAL"
               "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL."
               "For full copyright and/or warranty information see --copyright.";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, with_citation)
{
    argument_parser citation("test_parser", 2, argv1);
    citation.info.citation = "citation";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(citation.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_options_str +
               basic_version_str +
               "LEGAL"
               "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL."
               "In your academic works please cite: citation";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, empty_advanced_help)
{
    // Empty help call with -hh
    argument_parser parser2{"test_parser", 2, argv2};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser2.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_options_str +
               basic_version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, empty_version_call)
{
    // Empty version call
    argument_parser parser3{"test_parser", 2, argv3};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser3.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, version_call)
{
    // Version call with url and options.
    argument_parser parser4{"test_parser", 2, argv3};
    parser4.info.url = "www.seqan.de";
    parser4.add_option(option_value, 'i', "int", "this is a int option.");
    parser4.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser4.add_positional_option(pos_opt_value, "this is a positional option.");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser4.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_version_str +
               "URL"
               "www.seqan.de";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, do_not_print_hidden_options)
{
    // Add an option and request help.
    argument_parser parser5{"test_parser", 2, argv1};
    parser5.add_option(option_value, 'i', "int", "this is a int option.", option_spec::HIDDEN);
    parser5.add_flag(flag_value, 'f', "flag", "this is a flag.", option_spec::HIDDEN);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser5.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "===========" +
               basic_options_str +
               basic_version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, full_information)
{
    // Add synopsis, description, short description, positional option, option, flag, and example.
    argument_parser parser6{"test_parser", 2, argv1};
    parser6.info.synopsis.push_back("./some_binary_name synopsis");
    parser6.info.synopsis.push_back("./some_binary_name synopsis2");
    parser6.info.description.push_back("description");
    parser6.info.description.push_back("description2");
    parser6.info.short_description = "so short";
    parser6.add_option(option_value, 'i', "int", "this is a int option.");
    parser6.add_section("Flags");
    parser6.add_subsection("SubFlags");
    parser6.add_line("here come all the flags");
    parser6.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser6.add_positional_option(pos_opt_value, "this is a positional option.");
    parser6.info.examples.push_back("example");
    parser6.info.examples.push_back("example2");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser6.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser - so short\n"
               "======================\n"
               "SYNOPSIS\n"
               "./some_binary_name synopsis\n"
               "./some_binary_name synopsis2\n"
               "DESCRIPTION\n"
               "description\n"
               "description2\n"
               "POSITIONAL ARGUMENTS\n"
               "ARGUMENT-1 (List of std::string's)\n"
               "this is a positional option. Default: [].\n" +
               basic_options_str +
               "-i, --int (signed 32 bit integer)\n"
               "this is a int option. Default: 5.\n"
               "FLAGS\n"
               "SubFlags\n"
               "here come all the flags\n"
               "-f, --flag\n"
               "this is a flag.\n"
               "EXAMPLES\n"
               "example\n"
               "example2" +
               basic_version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, copyright)
{
    // Tests the --copyright call.
    const char * argvCopyright[] = {"./copyright", "--copyright"};
    argument_parser copyright("myApp", 2, argvCopyright);

    std::ifstream license_file{std::string{{SEQAN_INCLUDE_DIR}} + "/../LICENSE.md"};
    std::ranges::subrange<std::istreambuf_iterator<char>, std::istreambuf_iterator<char>> sub
    {
        std::istreambuf_iterator<char>(license_file),
        std::istreambuf_iterator<char>()
    };

    detail::consume(sub | view::take_until_and_consume(is_char<'`'>));
    std::string license_string{sub | view::drop(1) | view::take_until(is_char<'`'>)};

    // Test --copyright with empty short and long copyright info.
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(copyright.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std_cout = testing::internal::GetCapturedStdout();

        expected = "================================================================================\n"
                   "Copyright information for myApp:\n"
                   "--------------------------------------------------------------------------------\n"
                   "myApp copyright information not available.\n"
                   "================================================================================\n"
                   "This program contains SeqAn3 code licensed under the following terms:\n"
                   "--------------------------------------------------------------------------------\n"
                   + license_string;

        EXPECT_EQ(std_cout, expected);
    }

    // Test --copyright with a non-empty short copyright and an empty long copyright.
    copyright.info.short_copyright = "short copyright line 1\nshort copyright line 2";
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(copyright.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std_cout = testing::internal::GetCapturedStdout();

        expected = "================================================================================\n"
                   "Copyright information for myApp:\n"
                   "--------------------------------------------------------------------------------\n"
                   "myApp full copyright information not available. Displaying short copyright information instead:\n"
                   "short copyright line 1\n"
                   "short copyright line 2\n"
                   "================================================================================\n"
                   "This program contains SeqAn3 code licensed under the following terms:\n"
                   "--------------------------------------------------------------------------------\n"
                   + license_string;

        EXPECT_EQ(std_cout, expected);
    }

    // Test --copyright with a non-empty short copyright and a non-empty long copyright.
    copyright.info.long_copyright = "long copyright line 1\nlong copyright line 2";
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(copyright.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std_cout = testing::internal::GetCapturedStdout();

        expected = "================================================================================\n"
                   "Copyright information for myApp:\n"
                   "--------------------------------------------------------------------------------\n"
                   "long copyright line 1\n"
                   "long copyright line 2\n"
                   "================================================================================\n"
                   "This program contains SeqAn3 code licensed under the following terms:\n"
                   "--------------------------------------------------------------------------------\n"
                   + license_string;

        EXPECT_EQ(std_cout, expected);
    }
}
