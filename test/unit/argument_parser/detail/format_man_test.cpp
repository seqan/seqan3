// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_man.hpp>

using namespace seqan3;

// Reused global variables
int option_value{5};
bool flag_value{};
int8_t non_list_pos_opt_value{1};
std::vector<std::string> list_pos_opt_value{};
const char * argv[] = {"./format_man_test --version-check 0", "--export-help", "man"};
std::string my_stdout{};
std::string expected{".TH DEFAULT 1 \"December 01, 1994\" \"default 01.01.01\" \"default_man_page_title\"\n"
                     ".SH NAME\n"
                     "default \\- A short description here.\n"
                     ".SH SYNOPSIS\n"
                     "\\fB./format_man_test\\fP synopsis\n"
                     ".br\n"
                     "\\fB./format_man_test\\fP synopsis2\n"
                     ".SH DESCRIPTION\n"
                     "description\n"
                     ".sp\n"
                     "description2\n"
                     ".SH POSITIONAL ARGUMENTS\n"
                     ".TP\n"
                     "\\fBARGUMENT-1\\fP (\\fIsigned 8 bit integer\\fP)\n"
                     "this is a positional option. \n"
                     ".TP\n"
                     "\\fBARGUMENT-2\\fP (\\fIList\\fP of \\fIstd::string\\fP's)\n"
                     "this is a positional option. Default: []. \n"
                     ".SH OPTIONS\n"
                     ".SS Basic options:\n"
                     ".TP\n"
                     "\\fB-h\\fP, \\fB--help\\fP\n"
                     "Prints the help page.\n"
                     ".TP\n"
                     "\\fB-hh\\fP, \\fB--advanced-help\\fP\n"
                     "Prints the help page including advanced options.\n"
                     ".TP\n"
                     "\\fB--version\\fP\n"
                     "Prints the version information.\n"
                     ".TP\n"
                     "\\fB--copyright\\fP\n"
                     "Prints the copyright/license information.\n"
                     ".TP\n"
                     "\\fB--export-help\\fP (std::string)\n"
                     "Export the help page information. Value must be one of [html, man].\n"
                     ".TP\n"
                     "\\fB--version-check\\fP (bool)\n"
                     "Whether to to check for the newest app version. Default: 1.\n"
                     ".SS \n"
                     ".TP\n"
                     "\\fB-i\\fP, \\fB--int\\fP (\\fIsigned 32 bit integer\\fP)\n"
                     "this is a int option. Default: 5. \n"
                     ".TP\n"
                     "\\fB-j\\fP, \\fB--jint\\fP (\\fIsigned 32 bit integer\\fP)\n"
                     "this is a required int option. \n"
                     ".SH FLAGS\n"
                     ".SS SubFlags\n"
                     "here come all the flags\n"
                     ".TP\n"
                     "\\fB-f\\fP, \\fB--flag\\fP\n"
                     "this is a flag.\n"
                     ".TP\n"
                     "\\fB-k\\fP, \\fB--kflag\\fP\n"
                     "this is a flag.\n"
                     ".SH EXAMPLES\n"
                     "example\n"
                     ".sp\n"
                     "example2\n"};

// Full info parser initialisation
void dummy_init(argument_parser & parser)
{
    parser.info.date = "December 01, 1994";
    parser.info.version = "01.01.01";
    parser.info.man_page_title = "default_man_page_title";
    parser.info.short_description = "A short description here.";
    parser.info.synopsis.push_back("./format_man_test synopsis");
    parser.info.synopsis.push_back("./format_man_test synopsis2");
    parser.info.description.push_back("description");
    parser.info.description.push_back("description2");
    parser.add_option(option_value, 'i', "int", "this is a int option.");
    parser.add_option(option_value, 'j', "jint", "this is a required int option.", option_spec::REQUIRED);
    parser.add_section("Flags");
    parser.add_subsection("SubFlags");
    parser.add_line("here come all the flags");
    parser.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser.add_flag(flag_value, 'k', "kflag", "this is a flag.");
    parser.add_positional_option(non_list_pos_opt_value, "this is a positional option.");
    parser.add_positional_option(list_pos_opt_value, "this is a positional option.");
    parser.info.examples.push_back("example");
    parser.info.examples.push_back("example2");
}

TEST(format_man, empty_information)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};
    parser.info.date = "December 01, 1994";
    parser.info.version = "01.01.01";
    parser.info.man_page_title = "default_man_page_title";
    parser.info.short_description = "A short description here.";

    std::string expected_short{".TH DEFAULT 1 \"December 01, 1994\" \"default 01.01.01\" \"default_man_page_title\"\n"
                               ".SH NAME\n"
                               "default \\- A short description here.\n"
                               ".SH OPTIONS\n"
                               ".SS Basic options:\n"
                               ".TP\n"
                               "\\fB-h\\fP, \\fB--help\\fP\n"
                               "Prints the help page.\n"
                               ".TP\n"
                               "\\fB-hh\\fP, \\fB--advanced-help\\fP\n"
                               "Prints the help page including advanced options.\n"
                               ".TP\n"
                               "\\fB--version\\fP\n"
                               "Prints the version information.\n"
                               ".TP\n"
                               "\\fB--copyright\\fP\n"
                               "Prints the copyright/license information.\n"
                               ".TP\n"
                               "\\fB--export-help\\fP (std::string)\n"
                               "Export the help page information. Value must be one of [html, man].\n"
                               ".TP\n"
                               "\\fB--version-check\\fP (bool)\n"
                               "Whether to to check for the newest app version. Default: 1.\n"
                               ".SS \n"};

    // Test the dummy parser with minimal information.
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected_short);
}

TEST(format_man, full_information)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Test the dummy parser without any copyright or citations.
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST(format_man, full_info_short_copyright)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Add a short copyright and test the dummy parser.
    parser.info.short_copyright = "short copyright";
    expected += std::string{".SH LEGAL\n"
                            "\\fBdefault Copyright:\\fR short copyright\n"
                            ".br\n"
                            "\\fBSeqAn Copyright:\\fR 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.\n"
                            ".br\n"};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST(format_man, full_info_short_and_citation)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Add a short copyright & citation and test the dummy parser.
    parser.info.short_copyright = "short copyright";
    parser.info.citation = "citation";
    expected += std::string{"\\fBIn your academic works please cite:\\fR citation\n"
                            ".br\n"};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST(format_man, full_info_short_long_and_citation)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Add a short copyright & citation & long copyright and test the dummy parser.
    parser.info.short_copyright = "short copyright";
    parser.info.citation = "citation";
    parser.info.long_copyright = "looong copyright";
    expected += std::string{"For full copyright and/or warranty information see \\fB--copyright\\fR.\n"};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}
