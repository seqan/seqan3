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

const char * argv[] = {"./ctd_add_test", 
                       "--export-help", 
                       "ctd"};

TEST(ctd_format, empty_information)
{
    // Test '--export-help=ctd' call.
    argument_parser parser0{"empty_options", 
                            3, 
                            argv};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser0.parse(), 
                ::testing::ExitedWithCode(EXIT_SUCCESS), 
                "");
    EXPECT_EQ(testing::internal::GetCapturedStdout(),
              "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
              "<tool name=\"empty_options\" version=\"0.0.0.0\" ctdVersion=\"1.7.0\">\n"
              "\t<description/>\n"
              "\t<manual/>\n"
              "\t<cli/>\n"
              "\t<PARAMETERS version=\"1.7.0\">\n"
              "\t\t<NODE name=\"empty_options\" description=\"\"/>\n"
              "\t</PARAMETERS>\n"
              "</tool>\n"
              "\n");
}

TEST (ctd_format, test_valid_app_name)
{
    // App name cannot contain space characters.
    argument_parser parser0{"empty options",
                            3,
                            argv};
    EXPECT_THROW(parser0.parse(), 
                 parser_design_error);
    
    // App name cannot contain non-alphanumeric characters different from -/_
    argument_parser parser1{"empty.options",
                            3,
                            argv};
    EXPECT_THROW(parser1.parse(), 
                 parser_design_error);
}

TEST (ctd_format, test_add_option)
{
    std::string opt_a{};
    std::string opt_b{};
    std::string opt_c{};
    std::string opt_d{};

    argument_parser parser{"test_add_option",
                           3,
                           argv}; 
    parser.add_option(opt_a, 
                      'a', 
                      "", 
                      "Description short option A");
    parser.add_option(opt_b, 
                      'b', 
                      "option-b", 
                      "Description long option B",
                      ADVANCED);
    parser.add_option(opt_c,
                      'c',
                      "",
                      "Description short option C",
                      REQUIRED);
    parser.add_option(opt_d,
                      'd',
                      "",
                      "Description short option D",
                      HIDDEN);
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), 
                ::testing::ExitedWithCode(EXIT_SUCCESS), 
                "");
    EXPECT_EQ(testing::internal::GetCapturedStdout(),
              "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
              "<tool name=\"test_add_option\" version=\"0.0.0.0\" ctdVersion=\"1.7.0\">\n"
              "\t<description/>\n"
              "\t<manual/>\n"
              "\t<cli>\n"
              "\t\t<clielement optionIdentifier=\"-a\" isList=\"false\">\n"
              "\t\t\t<mapping referenceName=\"test_add_option.a\"/>\n"
              "\t\t</clielement>\n"
              "\t\t<clielement optionIdentifier=\"--option-b\" isList=\"false\">\n"
              "\t\t\t<mapping referenceName=\"test_add_option.option-b\"/>\n"
              "\t\t</clielement>\n"
              "\t\t<clielement optionIdentifier=\"-c\" isList=\"false\">\n"
              "\t\t\t<mapping referenceName=\"test_add_option.c\"/>\n"
              "\t\t</clielement>\n"
              "\t</cli>\n"
              "\t<PARAMETERS version=\"1.7.0\">\n"
              "\t\t<NODE name=\"test_add_option\" description=\"\">\n"
              "\t\t\t<ITEM name=\"a\" type=\"string\" description=\"Description short option A\" restrictions=\"\" required=\"false\" advanced=\"false\" value=\"\"/>\n"
              "\t\t\t<ITEM name=\"option-b\" type=\"string\" description=\"Description long option B\" restrictions=\"\" required=\"false\" advanced=\"true\" value=\"\"/>\n"
              "\t\t\t<ITEM name=\"c\" type=\"string\" description=\"Description short option C\" restrictions=\"\" required=\"true\" advanced=\"false\" value=\"\"/>\n"
              "\t\t</NODE>\n"
              "\t</PARAMETERS>\n"
              "</tool>\n"
              "\n");
}

TEST (ctd_format, test_add_flag) 
{
    bool opt_a{};

    argument_parser parser{"test_add_flag",
                           3,
                           argv}; 
    parser.add_flag(opt_a, 
                    'a', 
                    "", 
                    "Description short option A");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), 
                ::testing::ExitedWithCode(EXIT_SUCCESS), 
                "");
    EXPECT_EQ(testing::internal::GetCapturedStdout(),
              "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
              "<tool name=\"test_add_flag\" version=\"0.0.0.0\" ctdVersion=\"1.7.0\">\n"
              "\t<description/>\n"
              "\t<manual/>\n"
              "\t<cli>\n"
              "\t\t<clielement optionIdentifier=\"-a\" isList=\"false\">\n"
              "\t\t\t<mapping referenceName=\"test_add_flag.a\"/>\n"
              "\t\t</clielement>\n"
              "\t</cli>\n"
              "\t<PARAMETERS version=\"1.7.0\">\n"
              "\t\t<NODE name=\"test_add_flag\" description=\"\">\n"
              "\t\t\t<ITEM name=\"a\" type=\"bool\" description=\"Description short option A\" restrictions=\"\" required=\"false\" advanced=\"false\" value=\"\"/>\n"
              "\t\t</NODE>\n"
              "\t</PARAMETERS>\n"
              "</tool>\n"
              "\n");
}

TEST (ctd_format, test_add_positional_option)
{
    std::string arg_a{};
    std::string arg_b{};

    argument_parser parser{"test_add_positional_option",
                           3,
                           argv}; 
    parser.add_positional_option(arg_a, 
                                 "Description argument A");
    parser.add_positional_option(arg_b, 
                                 "Description argument B");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser.parse(), 
                ::testing::ExitedWithCode(EXIT_SUCCESS), 
                "");
    EXPECT_EQ(testing::internal::GetCapturedStdout(),
              "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
              "<tool name=\"test_add_positional_option\" version=\"0.0.0.0\" ctdVersion=\"1.7.0\">\n"
              "\t<description/>\n"
              "\t<manual/>\n"
              "\t<cli>\n"
              "\t\t<clielement optionIdentifier=\"\" isList=\"false\">\n"
              "\t\t\t<mapping referenceName=\"test_add_positional_option.argument-0\"/>\n"
              "\t\t</clielement>\n"
              "\t\t<clielement optionIdentifier=\"\" isList=\"false\">\n"
              "\t\t\t<mapping referenceName=\"test_add_positional_option.argument-1\"/>\n"
              "\t\t</clielement>\n"
              "\t</cli>\n"
              "\t<PARAMETERS version=\"1.7.0\">\n"
              "\t\t<NODE name=\"test_add_positional_option\" description=\"\">\n"
              "\t\t\t<ITEM name=\"argument-0\" type=\"string\" description=\"Description argument A\" restrictions=\"\" required=\"true\" advanced=\"false\" value=\"\"/>\n"
              "\t\t\t<ITEM name=\"argument-1\" type=\"string\" description=\"Description argument B\" restrictions=\"\" required=\"true\" advanced=\"false\" value=\"\"/>\n"
              "\t\t</NODE>\n"
              "\t</PARAMETERS>\n"
              "</tool>\n"
              "\n");
}
