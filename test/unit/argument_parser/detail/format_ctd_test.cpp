// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_ctd.hpp>

using namespace seqan3;

const char * argv0[] = {"./ctd_add_test", 
                        "--export-help", 
                        "ctd"};
const char * argv1[] = {"./ctd_add_test", 
                        "--export-help=ctd"};

TEST(ctd_format, empty_information)
{
    // Test the content of a CTD file when no command line option are configured.
    std::string expected_ctd_file = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                                    "<tool name=\"empty_options\" version=\"0.0.0.0\" ctdVersion=\"1.7.0\">\n"
                                    "\t<description/>\n"
                                    "\t<manual/>\n"
                                    "\t<cli/>\n"
                                    "\t<PARAMETERS version=\"1.7.0\">\n"
                                    "\t\t<NODE name=\"empty_options\" description=\"\"/>\n"
                                    "\t</PARAMETERS>\n"
                                    "</tool>\n"
                                    "\n";

    // Test '--export-help=ctd' call.
    argument_parser parser0{"empty_options", 
                            3, 
                            argv0};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser0.parse(), 
                ::testing::ExitedWithCode(EXIT_SUCCESS), 
                "");
    EXPECT_EQ(testing::internal::GetCapturedStdout(),
              expected_ctd_file);
    
    // Test '--export-help ctd' call.
    argument_parser parser1{"empty_options", 
                            2, 
                            argv1};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), 
                ::testing::ExitedWithCode(EXIT_SUCCESS), 
                "");
    EXPECT_EQ(testing::internal::GetCapturedStdout(),
              expected_ctd_file);
}

TEST (ctd_test, test_valid_app_name)
{
    // App name cannot contain space characters.
    argument_parser parser0{"empty options",
                            3,
                            argv0};
    EXPECT_THROW(parser0.parse(), 
                 parser_design_error);
    
    // App name cannot contain non-alphanumeric characters different from -/_
    argument_parser parser1{"empty.options",
                            3,
                            argv0};
    EXPECT_THROW(parser1.parse(), 
                 parser_design_error);
}

