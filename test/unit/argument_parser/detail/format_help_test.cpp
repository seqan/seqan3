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
    EXPECT_THROW(parser0.parse(), parser_interruption);
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
    EXPECT_THROW(parser1.parse(), parser_interruption);
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
    EXPECT_THROW(short_copy.parse(), parser_interruption);
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
    EXPECT_THROW(long_copy.parse(), parser_interruption);
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
    EXPECT_THROW(citation.parse(), parser_interruption);
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
    // EXPECT_THROW(copyright.parse(), parser_interruption);
    // stdout = testing::internal::GetCapturedStdout();
    // expected = "whatever is expected";
    // EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               // expected | ranges::view::remove_if(is_space)));

    // Empty help call with -hh
    const char * argv2[] = {"./help_add_test", "-hh"};
    argument_parser parser2("test_parser_2", 2, argv2);
    testing::internal::CaptureStdout();
    EXPECT_THROW(parser2.parse(), parser_interruption);
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
    EXPECT_THROW(parser3.parse(), parser_interruption);
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
    EXPECT_THROW(parser4.parse(), parser_interruption);
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
    EXPECT_THROW(parser5.parse(), parser_interruption);
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
    EXPECT_THROW(parser6.parse(), parser_interruption);
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

   // EXPECT_THROW(parser6.parse(), parser_interruption);
}
