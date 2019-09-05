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
struct format_man_test : public ::testing::Test
{
    int option_value{5};
    bool flag_value{};
    int8_t non_list_pos_opt_value{1};
    std::vector<std::string> list_pos_opt_value{};
    std::string my_stdout{};
    const char * argv[4] = {"./format_man_test --version-check 0", "--export-help", "man"};
    std::string expected = R"(.TH DEFAULT 1 "December 01, 1994" "default 01.01.01" "default_man_page_title"
.SH NAME
default \- A short description here.
.SH SYNOPSIS
\fB./format_man_test\fP synopsis
.br
\fB./format_man_test\fP synopsis2
.SH DESCRIPTION
description
.sp
description2
.SH POSITIONAL ARGUMENTS
.TP
\fBARGUMENT-1\fP (\fIsigned 8 bit integer\fP)
this is a positional option. 
.TP
\fBARGUMENT-2\fP (\fIList\fP of \fIstd::string\fP's)
this is a positional option. Default: []. 
.SH OPTIONS
.SS Basic options:
.TP
\fB-h\fP, \fB--help\fP
Prints the help page.
.TP
\fB-hh\fP, \fB--advanced-help\fP
Prints the help page including advanced options.
.TP
\fB--version\fP
Prints the version information.
.TP
\fB--copyright\fP
Prints the copyright/license information.
.TP
\fB--export-help\fP (std::string)
Export the help page information. Value must be one of [html, man].
.TP
\fB--version-check\fP (bool)
Whether to to check for the newest app version. Default: 1.
.SS 
.TP
\fB-i\fP, \fB--int\fP (\fIsigned 32 bit integer\fP)
this is a int option. Default: 5. 
.TP
\fB-j\fP, \fB--jint\fP (\fIsigned 32 bit integer\fP)
this is a required int option. 
.SH FLAGS
.SS SubFlags
here come all the flags
.TP
\fB-f\fP, \fB--flag\fP
this is a flag.
.TP
\fB-k\fP, \fB--kflag\fP
this is a flag.
.SH EXAMPLES
example
.sp
example2
)";

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
};

TEST_F(format_man_test, empty_information)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};
    parser.info.date = "December 01, 1994";
    parser.info.version = "01.01.01";
    parser.info.man_page_title = "default_man_page_title";
    parser.info.short_description = "A short description here.";

    std::string expected_short = R"(.TH DEFAULT 1 "December 01, 1994" "default 01.01.01" "default_man_page_title"
.SH NAME
default \- A short description here.
.SH OPTIONS
.SS Basic options:
.TP
\fB-h\fP, \fB--help\fP
Prints the help page.
.TP
\fB-hh\fP, \fB--advanced-help\fP
Prints the help page including advanced options.
.TP
\fB--version\fP
Prints the version information.
.TP
\fB--copyright\fP
Prints the copyright/license information.
.TP
\fB--export-help\fP (std::string)
Export the help page information. Value must be one of [html, man].
.TP
\fB--version-check\fP (bool)
Whether to to check for the newest app version. Default: 1.
.SS 
)";

    // Test the dummy parser with minimal information.
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected_short);
}

TEST_F(format_man_test, full_information)
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

TEST_F(format_man_test, full_info_short_copyright)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Add a short copyright and test the dummy parser.
    parser.info.short_copyright = "short copyright";
    expected += R"(.SH LEGAL
\fBdefault Copyright:\fR short copyright
.br
\fBSeqAn Copyright:\fR 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.
.br
)";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST_F(format_man_test, full_info_short_and_citation)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Add a short copyright & citation and test the dummy parser.
    parser.info.short_copyright = "short copyright";
    parser.info.citation = "citation";
    expected += R"(.SH LEGAL
\fBdefault Copyright:\fR short copyright
.br
\fBSeqAn Copyright:\fR 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.
.br
\fBIn your academic works please cite:\fR citation
.br
)";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST_F(format_man_test, full_info_short_long_and_citation)
{
    // Create the dummy parser.
    argument_parser parser{"default", 3, argv};

    // Fill out the dummy parser with options and flags and sections and subsections.
    dummy_init(parser);
    // Add a short copyright & citation & long copyright and test the dummy parser.
    parser.info.short_copyright = "short copyright";
    parser.info.citation = "citation";
    parser.info.long_copyright = "looong copyright";
    expected += R"(.SH LEGAL
\fBdefault Copyright:\fR short copyright
.br
\fBSeqAn Copyright:\fR 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.
.br
\fBIn your academic works please cite:\fR citation
.br
For full copyright and/or warranty information see \fB--copyright\fR.
)";
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}
