// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/remove_if.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_help.hpp>
#include <seqan3/io/stream/parse_condition.hpp>

using namespace seqan3;

TEST(help_add_test, add_option)
{
    std::string stdout;
    std::string expected;
    int option_value;
    bool flag_value;
    std::vector<std::string> pos_opt_value;

    // Empty call with no options given. For detail::format_short_help
    const char * argv0[] = {"./help_add_test"};
    argument_parser parser0("empty_options", 1, argv0);
    parser0.info.synopsis.push_back("synopsis");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser0.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("empty_options"
                           "============="
                           "empty_options synopsis"
                           "Try -h or --help for more information.");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Empty help call with -h
    const char * argv1[] = {"./help_add_test", "-h"};
    argument_parser parser1("test_parser", 2, argv1);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("test_parser"
                                       "==========="
                                       "VERSION"
                                       "Last update:"
                                       "test_parser version:"
                                       "SeqAn version: 3.0.0");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Again, but with short copyright, long copyright, and citation.
    argument_parser short_copy("test_parser", 2, argv1);
    short_copy.info.short_copyright = "short";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(short_copy.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("test_parser"
                           "==========="
                           "VERSION"
                           "Last update:"
                           "test_parser version:"
                           "SeqAn version: 3.0.0"
                           "LEGAL"
                           "test_parser Copyright: short"
                           "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    argument_parser long_copy("test_parser", 2, argv1);
    long_copy.info.long_copyright = "long";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(long_copy.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("test_parser"
                           "==========="
                           "VERSION"
                           "Last update:"
                           "test_parser version:"
                           "SeqAn version: 3.0.0"
                           "LEGAL"
                           "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL."
                           "For full copyright and/or warranty information see --copyright.");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    argument_parser citation("test_parser", 2, argv1);
    citation.info.citation = "citation";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(citation.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("test_parser"
                           "==========="
                           "VERSION"
                           "Last update:"
                           "test_parser version:"
                           "SeqAn version: 3.0.0"
                           "LEGAL"
                           "SeqAn Copyright: 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL."
                           "In your academic works please cite: citation");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Tests the --copyright call. TODO: should be uncommented/fixed after --copyright is actually implemented.
    // const char * argvCopyright[] = {"./copyright", "--copyright"};
    // argument_parser copyright("copyright", 2, argvCopyright);
    // copyright.info.long_copyright = "long copyright";
    // testing::internal::CaptureStdout();
    // EXPECT_EXIT(copyright.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    // stdout = testing::internal::GetCapturedStdout();
    // expected = "whatever is expected";
    // EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               // expected | ranges::view::remove_if(is_space)));

    // Empty help call with -hh
    const char * argv2[] = {"./help_add_test", "-hh"};
    argument_parser parser2("test_parser_2", 2, argv2);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser2.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("test_parser_2"
                           "============="
                           "VERSION"
                           "Last update:"
                           "test_parser_2 version:"
                           "SeqAn version: 3.0.0");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Empty version call
    const char * argv3[] = {"./help_add_test", "--version"};
    argument_parser parser3("version", 2, argv3);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser3.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("version"
                           "======="
                           "VERSION"
                           "Last update:"
                           "version version:"
                           "SeqAn version: 3.0.0");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Version call with url and options.
    argument_parser parser4("versionURL", 2, argv3);
    parser4.info.url = "www.seqan.de";
    parser4.add_option(option_value, 'i', "int", "this is a int option.");
    parser4.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser4.add_positional_option(pos_opt_value, "this is a positional option.");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser4.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("versionURL"
                           "=========="
                           "VERSION"
                           "Last update:"
                           "versionURL version:"
                           "SeqAn version: 3.0.0"
                           "URL"
                           "www.seqan.de");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Add an option and request help.
    argument_parser parser5("hidden", 2, argv1);
    parser5.add_option(option_value, 'i', "int", "this is a int option.", option_spec::HIDDEN);
    parser5.add_flag(flag_value, 'f', "flag", "this is a flag.", option_spec::HIDDEN);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser5.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("hidden"
                           "======"
                           "OPTIONS"
                           "VERSION"
                           "Last update:"
                           "hidden version:"
                           "SeqAn version: 3.0.0");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

    // Add synopsis, description, short description, positional option, option, flag, and example.
    argument_parser parser6("full", 2, argv1);
    parser6.info.synopsis.push_back("synopsis");
    parser6.info.synopsis.push_back("synopsis2");
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
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("full - so short"
                           "==============="
                           "SYNOPSIS"
                           "full synopsis"
                           "full synopsis2"
                           "DESCRIPTION"
                           "description"
                           "description2"
                           "POSITIONAL ARGUMENTS"
                           "ARGUMENT-1 List of STRING's"
                           "this is a positional option."
                           "OPTIONS"
                           "-i, --int INT (32 bit)"
                           "this is a int option."
                           "-f, --flag"
                           "this is a flag."
                           "EXAMPLES"
                           "example"
                           "example2"
                           "VERSION"
                           "Last update:"
                           "full version:"
                           "SeqAn version: 3.0.0");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

   // EXPECT_EXIT(parser6.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
}
