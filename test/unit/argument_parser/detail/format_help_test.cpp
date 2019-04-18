// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <fstream>

#include <gtest/gtest.h>

#include <seqan3/std/ranges>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_help.hpp>
#include <seqan3/io/stream/parse_condition.hpp>

using namespace seqan3;

// reused global variables
std::string std_cout;
std::string expected;
int option_value;
bool flag_value;
std::vector<std::string> pos_opt_value;
const char * argv0[] = {"./help_add_test"};
const char * argv1[] = {"./help_add_test", "-h"};
const char * argv2[] = {"./help_add_test", "-hh"};
const char * argv3[] = {"./help_add_test", "--version"};

std::string version_str = std::to_string(SEQAN3_VERSION_MAJOR) + "."
                          + std::to_string(SEQAN3_VERSION_MINOR) + "."
                          + std::to_string(SEQAN3_VERSION_PATCH);

TEST(help_page_printing, short_help)
{
    // Empty call with no options given. For detail::format_short_help
    argument_parser parser0("empty_options", 1, argv0);
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
    argument_parser parser1("test_parser", 2, argv1);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "==========="
               "VERSION"
               "Last update:"
               "test_parser version:"
               "SeqAn version: " + version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, with_copyright)
{
    // Again, but with copyright, license, and citation.
    argument_parser copyright("test_parser", 2, argv1);
    copyright.info.copyright = "copyright";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(copyright.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "==========="
               "VERSION"
               "Last update:"
               "test_parser version:"
               "SeqAn version: " + version_str +
               "LEGAL"
               "test_parser Copyright: copyright"
               "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, with_license)
{
    argument_parser license("test_parser", 2, argv1);
    license.info.license = "license";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(license.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser"
               "==========="
               "VERSION"
               "Last update:"
               "test_parser version:"
               "SeqAn version: " + version_str +
               "LEGAL"
               "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL."
               "For full license and/or warranty information see --license.";
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
               "==========="
               "VERSION"
               "Last update:"
               "test_parser version:"
               "SeqAn version: " + version_str +
               "LEGAL"
               "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL."
               "In your academic works please cite: citation";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, empty_advanced_help)
{
    // Empty help call with -hh
    argument_parser parser2("test_parser_2", 2, argv2);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser2.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "test_parser_2"
               "============="
               "VERSION"
               "Last update:"
               "test_parser_2 version:"
               "SeqAn version: " + version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, empty_version_call)
{
    // Empty version call
    argument_parser parser3("version", 2, argv3);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser3.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "version"
               "======="
               "VERSION"
               "Last update:"
               "version version:"
               "SeqAn version: " + version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, version_call)
{
    // Version call with url and options.
    argument_parser parser4("versionURL", 2, argv3);
    parser4.info.url = "www.seqan.de";
    parser4.add_option(option_value, 'i', "int", "this is a int option.");
    parser4.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser4.add_positional_option(pos_opt_value, "this is a positional option.");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser4.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "versionURL"
               "=========="
               "VERSION"
               "Last update:"
               "versionURL version:"
               "SeqAn version: " + version_str +
               "URL"
               "www.seqan.de";
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, do_not_print_hidden_options)
{
    // Add an option and request help.
    argument_parser parser5("hidden", 2, argv1);
    parser5.add_option(option_value, 'i', "int", "this is a int option.", option_spec::HIDDEN);
    parser5.add_flag(flag_value, 'f', "flag", "this is a flag.", option_spec::HIDDEN);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser5.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "hidden"
               "======"
               "OPTIONS"
               "VERSION"
               "Last update:"
               "hidden version:"
               "SeqAn version: " + version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, full_information)
{
    // Add synopsis, description, short description, positional option, option, flag, and example.
    argument_parser parser6("full", 2, argv1);
    parser6.info.synopsis.push_back("./some_binary_name synopsis");
    parser6.info.synopsis.push_back("./some_binary_name synopsis2");
    parser6.info.description.push_back("description");
    parser6.info.description.push_back("description2");
    parser6.info.short_description = "so short";
    parser6.add_option(option_value, 'i', "int", "this is a int option.");
    parser6.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser6.add_positional_option(pos_opt_value, "this is a positional option.");
    parser6.info.examples.push_back("example");
    parser6.info.examples.push_back("example2");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser6.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std_cout = testing::internal::GetCapturedStdout();
    expected = "full - so short"
               "==============="
               "SYNOPSIS"
               "./some_binary_name synopsis"
               "./some_binary_name synopsis2"
               "DESCRIPTION"
               "description"
               "description2"
               "POSITIONAL ARGUMENTS"
               "ARGUMENT-1 (List of std::string's)"
               "this is a positional option."
               "OPTIONS"
               "-i, --int (signed 32 bit integer)"
               "this is a int option."
               "-f, --flag"
               "this is a flag."
               "EXAMPLES"
               "example"
               "example2"
               "VERSION"
               "Last update:"
               "full version:"
               "SeqAn version: " + version_str;
    EXPECT_TRUE(ranges::equal((std_cout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
}

TEST(help_page_printing, license)
{
    // Tests the --license call.
    const char * argvLicense[] = {"./copyright", "--license"};
    argument_parser license("myApp", 2, argvLicense);

    std::ifstream license_file{std::string{{SEQAN_INCLUDE_DIR}} + "/../LICENSE"};
    std::string license_string{(std::istreambuf_iterator<char>(license_file)),
                                std::istreambuf_iterator<char>()};

    // Test --license with empty copyright and license info.
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(license.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std_cout = testing::internal::GetCapturedStdout();

        expected = "================================================================================\n"
                   "License information for myApp:\n"
                   "--------------------------------------------------------------------------------\n"
                   "myApp copyright and license information not available.\n"
                   "================================================================================\n"
                   "This program contains SeqAn3 code licensed under the following terms:\n"
                   "--------------------------------------------------------------------------------\n"
                   + license_string;

        EXPECT_EQ(std_cout, expected);
    }

    // Test --copyright with a non-empty copyright and an empty license.
    license.info.copyright = "copyright line 1\ncopyright line 2";
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(license.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std_cout = testing::internal::GetCapturedStdout();

        expected = "================================================================================\n"
                   "License information for myApp:\n"
                   "--------------------------------------------------------------------------------\n"
                   "myApp full license information not available. Displaying copyright information instead:\n"
                   "copyright line 1\n"
                   "copyright line 2\n"
                   "================================================================================\n"
                   "This program contains SeqAn3 code licensed under the following terms:\n"
                   "--------------------------------------------------------------------------------\n"
                   + license_string;

        EXPECT_EQ(std_cout, expected);
    }

    // Test --license with a non-empty copyright and a non-empty license.
    license.info.license = "license line 1\nlicense line 2";
    {
        testing::internal::CaptureStdout();
        EXPECT_EXIT(license.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std_cout = testing::internal::GetCapturedStdout();

        expected = "================================================================================\n"
                   "License information for myApp:\n"
                   "--------------------------------------------------------------------------------\n"
                   "license line 1\n"
                   "license line 2\n"
                   "================================================================================\n"
                   "This program contains SeqAn3 code licensed under the following terms:\n"
                   "--------------------------------------------------------------------------------\n"
                   + license_string;

        EXPECT_EQ(std_cout, expected);
    }
}
