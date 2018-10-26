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
#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

TEST(validator_test, fullfill_concept)
{
    EXPECT_TRUE(validator_concept<detail::default_validator<int>>);
    EXPECT_TRUE(validator_concept<integral_range_validator<int>>);
    EXPECT_TRUE(validator_concept<value_list_validator<int>>);
    EXPECT_TRUE(validator_concept<regex_validator<std::string>>);
    EXPECT_TRUE(validator_concept<regex_validator<std::vector<std::string>>>);
    EXPECT_TRUE(validator_concept<detail::default_validator<std::vector<int>>>);
    EXPECT_TRUE(validator_concept<integral_range_validator<std::vector<int>>>);
    EXPECT_TRUE(validator_concept<value_list_validator<std::vector<int>>>);
    EXPECT_TRUE(validator_concept<file_ext_validator>);
}

TEST(validator_test, no_file)
{
    filesystem::path p{"./sandbox.fasta"};
    std::string s{"./stonebox.fasta"};
    file_existance_validator my_validator{};
    EXPECT_THROW(my_validator(p), parser_invalid_argument);
    EXPECT_THROW(my_validator(s), parser_invalid_argument);
}

TEST(validator_test, file_exists)
{
    test::tmp_filename tmp_file_name{""};
    std::ofstream tmp_file(tmp_file_name.get_path());
    file_existance_validator my_validator{};
    for (auto & file : filesystem::directory_iterator(tmp_file_name.get_path()))
        EXPECT_NO_THROW(my_validator(file));
    for (auto & file : filesystem::directory_iterator(tmp_file_name.get_path()))
          EXPECT_NO_THROW(my_validator(file));
}

TEST(validator_test, integral_range_validator_success)
{
    int option_value;
    std::vector<int> option_vector;

    // option
    const char * argv[] = {"./argument_parser_test", "-i", "10"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 'i', "int-option", "desc",
                      option_spec::DEFAULT, integral_range_validator<int>(1,20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // option - negative values
    const char * argv2[] = {"./argument_parser_test", "-i", "-10"};
    argument_parser parser2("test_parser", 3, argv2);
    parser2.add_option(option_value, 'i', "int-option", "desc",
                       option_spec::DEFAULT, integral_range_validator<int>(-20,20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "10"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value, "desc", integral_range_validator<int>(1,20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // positional option - negative values
    const char * argv4[] = {"./argument_parser_test", "--", "-10"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_value, "desc", integral_range_validator<int>(-20,20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    argument_parser parser5("test_parser", 5, argv5);
    parser5.add_option(option_vector, 'i', "int-option", "desc",
                       option_spec::DEFAULT, integral_range_validator<std::vector<int>>(-50,50));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 48);

    // positional option - vector
    option_vector.clear();
    const char * argv6[] = {"./argument_parser_test", "--", "-10", "1"};
    argument_parser parser6("test_parser", 4, argv6);
    parser6.add_positional_option(option_vector, "desc", integral_range_validator<std::vector<int>>(-20,20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser6.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 1);
}

TEST(validator_test, integral_range_validator_error)
{
    int option_value;
    std::vector<int> option_vector;

    // option - above max
    const char * argv[] = {"./argument_parser_test", "-i", "30"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 'i', "int-option", "desc",
                      option_spec::DEFAULT, integral_range_validator<int>(1,20));

    EXPECT_THROW(parser.parse(), validation_failed);

    // option - below min
    const char * argv2[] = {"./argument_parser_test", "-i", "-21"};
    argument_parser parser2("test_parser", 3, argv2);
    parser2.add_option(option_value, 'i', "int-option", "desc",
                       option_spec::DEFAULT, integral_range_validator<int>(-20,20));

    EXPECT_THROW(parser2.parse(), validation_failed);

    // positional option - above max
    const char * argv3[] = {"./argument_parser_test", "30"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value, "desc", integral_range_validator<int>(1,20));

    EXPECT_THROW(parser3.parse(), validation_failed);

    // positional option - below min
    const char * argv4[] = {"./argument_parser_test", "--", "-21"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_value, "desc", integral_range_validator<int>(-20,20));

    EXPECT_THROW(parser4.parse(), validation_failed);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-100"};
    argument_parser parser5("test_parser", 3, argv5);
    parser5.add_option(option_vector, 'i', "int-option", "desc",
                       option_spec::DEFAULT, integral_range_validator<std::vector<int>>(-50,50));

    EXPECT_THROW(parser5.parse(), validation_failed);

    // positional option - vector
    option_vector.clear();
    const char * argv6[] = {"./argument_parser_test", "--", "-10", "100"};
    argument_parser parser6("test_parser", 4, argv6);
    parser6.add_positional_option(option_vector, "desc", integral_range_validator<std::vector<int>>(-20,20));

    EXPECT_THROW(parser6.parse(), validation_failed);
}

TEST(validator_test, value_list_validator_success)
{
    std::string option_value;
    int option_value_int;
    std::vector<std::string> option_vector;
    std::vector<int> option_vector_int;

    // option
    const char * argv[] = {"./argument_parser_test", "-s", "ba"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 's', "string-option", "desc",
                      option_spec::DEFAULT, value_list_validator<std::string>({"ha", "ba", "ma"}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ba");

    // option with integers
    const char * argv2[] = {"./argument_parser_test", "-i", "-21"};
    argument_parser parser2("test_parser", 3, argv2);
    parser2.add_option(option_value_int, 'i', "int-option", "desc",
                       option_spec::DEFAULT, value_list_validator<int>({0, -21, 10}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value_int, -21);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "ma"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value, "desc", value_list_validator<std::string>({"ha", "ba", "ma"}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ma");

    // positional option - vector
    const char * argv4[] = {"./argument_parser_test", "ha", "ma"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_vector, "desc",
                                  value_list_validator<std::vector<std::string>>({"ha", "ba", "ma"}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "ha");
    EXPECT_EQ(option_vector[1], "ma");

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    argument_parser parser5("test_parser", 5, argv5);
    parser5.add_option(option_vector_int, 'i', "int-option", "desc",
                       option_spec::DEFAULT, value_list_validator<std::vector<int>>({-10,48,50}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector_int[0], -10);
    EXPECT_EQ(option_vector_int[1], 48);
}

TEST(validator_test, value_list_validator_error)
{
    std::string option_value;
    int option_value_int;
    std::vector<std::string> option_vector;
    std::vector<int> option_vector_int;

    // option
    const char * argv[] = {"./argument_parser_test", "-s", "sa"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 's', "string-option", "desc",
                      option_spec::DEFAULT, value_list_validator<std::string>({"ha", "ba", "ma"}));

    EXPECT_THROW(parser.parse(), validation_failed);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "30"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value_int, "desc", value_list_validator<int>({0, 5, 10}));

    EXPECT_THROW(parser3.parse(), validation_failed);

    // positional option - vector
    const char * argv4[] = {"./argument_parser_test", "fo", "ma"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_vector, "desc",
                                  value_list_validator<std::vector<std::string>>({"ha", "ba", "ma"}));

    EXPECT_THROW(parser4.parse(), validation_failed);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "488"};
    argument_parser parser5("test_parser", 5, argv5);
    parser5.add_option(option_vector_int, 'i', "int-option", "desc",
                       option_spec::DEFAULT, value_list_validator<std::vector<int>>({-10,48,50}));

    EXPECT_THROW(parser5.parse(), validation_failed);
}

TEST(validator_test, regex_validator_success)
{
    std::string option_value;
    std::vector<std::string> option_vector;
    regex_validator<std::string> email_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");
    regex_validator<std::vector<std::string>> email_vector_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");

    // option
    const char * argv[] = {"./argument_parser_test", "-s", "ballo@rollo.com"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 's', "string-option", "desc",
                      option_spec::DEFAULT, email_validator);

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ballo@rollo.com");

    // positional option
    const char * argv2[] = {"./argument_parser_test", "chr1"};
    argument_parser parser2("test_parser", 2, argv2);
    parser2.add_positional_option(option_value, "desc",
                                  regex_validator<std::string>{"^chr[0-9]+"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "chr1");

    // positional option - vector
    const char * argv3[] = {"./argument_parser_test", "rollo", "bollo", "lollo"};
    argument_parser parser3("test_parser", 4, argv3);
    parser3.add_positional_option(option_vector, "desc",
                                  regex_validator<std::vector<std::string>>{".*oll.*"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "rollo");
    EXPECT_EQ(option_vector[1], "bollo");
    EXPECT_EQ(option_vector[2], "lollo");

    // option - vector
    option_vector.clear();
    const char * argv4[] = {"./argument_parser_test", "-s", "rita@rambo.com", "-s", "tina@rambo.com"};
    argument_parser parser4("test_parser", 5, argv4);
    parser4.add_option(option_vector, 's', "string-option", "desc",
                       option_spec::DEFAULT, email_vector_validator);

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "rita@rambo.com");
    EXPECT_EQ(option_vector[1], "tina@rambo.com");

}

TEST(validator_test, regex_validator_error)
{
    std::string option_value;
    std::vector<std::string> option_vector;

    // option
    const char * argv[] = {"./argument_parser_test", "-s", "sally"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 's', "string-option", "desc",
                      option_spec::DEFAULT, regex_validator<std::string>("tt"));

    EXPECT_THROW(parser.parse(), validation_failed);

    // positional option
    const char * argv2[] = {"./argument_parser_test", "jessy"};
    argument_parser parser2("test_parser", 2, argv2);
    parser2.add_positional_option(option_value, "desc",
                                  regex_validator<std::string>("[0-9]"));

    EXPECT_THROW(parser2.parse(), validation_failed);

    // positional option - vector
    const char * argv3[] = {"./argument_parser_test", "rollo", "bttllo", "lollo"};
    argument_parser parser3("test_parser", 4, argv3);
    parser3.add_positional_option(option_vector, "desc",
                                  regex_validator<std::vector<std::string>>{".*oll.*"});

    EXPECT_THROW(parser3.parse(), validation_failed);

    // option - vector
    option_vector.clear();
    const char * argv4[] = {"./argument_parser_test", "-s", "gh", "-s", "tt"};
    argument_parser parser4("test_parser", 5, argv4);
    parser4.add_option(option_vector, 's', "string-option", "desc",
                       option_spec::DEFAULT, regex_validator<std::vector<std::string>>("tt"));

    EXPECT_THROW(parser4.parse(), validation_failed);
}
