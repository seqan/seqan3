// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <gtest/gtest.h>

#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

TEST(parse_test, parser_design_error)
{
    int option_value;

    // short option
    const char * argv[] = {"./argument_parser_test"};
    argument_parser parser("test_parser", 1, argv);
    parser.add_option(option_value, 'i', "int", "this is a int option.");
    EXPECT_THROW(parser.add_option(option_value, 'i', "aint", "oh oh same id."),
                 parser_design_error);

    // long option
    argument_parser parser2("test_parser", 1, argv);
    parser2.add_option(option_value, 'i', "int", "this is an int option.");
    EXPECT_THROW(parser2.add_option(option_value, 'a', "int", "oh oh another id."),
                 parser_design_error);

    // empty identifier
    argument_parser parser3("test_parser", 1, argv);
    EXPECT_THROW(parser3.add_option(option_value, '\0', "", "oh oh all is empty."),
                 parser_design_error);

    bool flag_value;

    // short flag
    argument_parser parser4("test_parser", 1, argv);
    parser4.add_flag(flag_value, 'i', "int1", "this is an int option.");
    EXPECT_THROW(parser4.add_flag(flag_value, 'i', "int2", "oh oh another id."),
                 parser_design_error);

    // long flag
    argument_parser parser5("test_parser", 1, argv);
    parser5.add_flag(flag_value, 'i', "int", "this is an int option.");
    EXPECT_THROW(parser5.add_flag(flag_value, 'a', "int", "oh oh another id."),
                 parser_design_error);

    // empty identifier
    argument_parser parser6("test_parser", 1, argv);
    EXPECT_THROW(parser6.add_flag(flag_value, '\0', "", "oh oh another id."),
                 parser_design_error);

    // positional option not at the end
    const char * argv2[] = {"./argument_parser_test", "arg1", "arg2", "arg3"};
    std::vector<int> vec;
    argument_parser parser7("test_parser", 4, argv2);
    parser7.add_positional_option(vec, "oh oh list not at the end.");
    parser7.add_positional_option(option_value, "desc.");
    EXPECT_THROW(parser7.parse(), parser_design_error);

    // using h, help, advanced-help, and export-help
    argument_parser parser8("test_parser", 1, argv);
    EXPECT_THROW(parser8.add_option(option_value, 'h', "", "-h is bad."),
                 parser_design_error);
    EXPECT_THROW(parser8.add_option(option_value, '\0', "help", "help is bad."),
                 parser_design_error);
    EXPECT_THROW(parser8.add_option(option_value, '\0', "advanced-help",
                 "advanced-help is bad"), parser_design_error);
    EXPECT_THROW(parser8.add_option(option_value, '\0', "export-help",
                 "export-help is bad"), parser_design_error);
}

TEST(parse_test, parse_called_twice)
{
    std::string option_value;

    const char * argv[] = {"./argument_parser_test", "-s", "option_string"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 's', "string-option", "this is a string option.");

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "option_string");

    EXPECT_THROW(parser.parse(), parser_design_error);
}
