// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>

#include <range/v3/view/remove_if.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

std::string const basic_options_str = "OPTIONS"
                                      "Basic options:"
                                      "-h, --help Prints the help page."
                                      "-hh, --advanced-help Prints the help page including advanced options."
                                      "--version Prints the version information."
                                      "--copyright Prints the copyright/license information."
                                      "--export-help (std::string) Export the help page information. "
                                                                   "Value must be one of [html, man].";

std::string const basic_version_str = "VERSION"
                                      "Last update:"
                                      "test_parser version:"
                                      "SeqAn version: " + seqan3_version;

TEST(validator_test, fullfill_concept)
{
    EXPECT_FALSE(validator_concept<int>);

    EXPECT_TRUE(validator_concept<detail::default_validator<int>>);
    EXPECT_TRUE(validator_concept<detail::default_validator<int> const>);
    EXPECT_TRUE(validator_concept<detail::default_validator<int> &>);

    EXPECT_TRUE(validator_concept<detail::default_validator<std::vector<int>>>);
    EXPECT_TRUE(validator_concept<arithmetic_range_validator>);
    EXPECT_TRUE(validator_concept<value_list_validator<double>>);
    EXPECT_TRUE(validator_concept<value_list_validator<std::string>>);
    EXPECT_TRUE(validator_concept<input_file_validator>);
    EXPECT_TRUE(validator_concept<output_file_validator>);
    EXPECT_TRUE(validator_concept<input_directory_validator>);
    EXPECT_TRUE(validator_concept<output_directory_validator>);
    EXPECT_TRUE(validator_concept<regex_validator>);

    EXPECT_TRUE(validator_concept<decltype(input_file_validator{{"t"}} | regex_validator{".*"})>);
}

TEST(validator_test, input_file)
{
    test::tmp_filename tmp_name{"testbox.fasta"};
    test::tmp_filename tmp_name_2{"testbox_2.fasta"};

    std::vector formats{std::string{"fa"}, std::string{"sam"}, std::string{"fasta"}};

    std::ofstream tmp_file(tmp_name.get_path());
    std::ofstream tmp_file_2(tmp_name_2.get_path());

    { // single file

        { // empty list of file.
            input_file_validator my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name.get_path()));
        }

        { // file already exists.
            std::filesystem::path does_not_exists{tmp_name.get_path()};
            does_not_exists.replace_extension(".bam");
            input_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exists), parser_invalid_argument);
        }

        { // file has wrong format.
            input_file_validator my_validator{std::vector{std::string{"sam"}}};
                EXPECT_THROW(my_validator(tmp_name.get_path()), parser_invalid_argument);
        }

        { // file has no extension.
            std::filesystem::path does_not_exists{tmp_name.get_path()};
            does_not_exists.replace_extension();
            input_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exists), parser_invalid_argument);
        }

        std::filesystem::path file_in_path;

        // option
        std::filesystem::path path = tmp_name.get_path();
        const char * argv[] = {"./argument_parser_test", "-i", path.c_str()};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(file_in_path, 'i', "int-option", "desc",
                          option_spec::DEFAULT, input_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(file_in_path.string(), path.string());
    }

    { // file list.
        std::vector<std::filesystem::path> input_files;

        // option
        std::filesystem::path path = tmp_name.get_path();
        std::filesystem::path path_2 = tmp_name_2.get_path();

        const char * argv[] = {"./argument_parser_test", path.c_str(), path_2.c_str()};
        argument_parser parser("test_parser", 3, argv);
        parser.add_positional_option(input_files, "desc", input_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(input_files.size(), 2u);
        EXPECT_EQ(input_files[0].string(), path.string());
        EXPECT_EQ(input_files[1].string(), path_2.string());
    }

    { // get help page message
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        argument_parser parser("test_parser", 2, argv);
        parser.add_positional_option(path, "desc", input_file_validator{formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser"
                               "==========="
                               "POSITIONAL ARGUMENTS"
                               "    ARGUMENT-1 (std::filesystem::path)"
                               "          desc Default: \"\". Input file formats: fa, sam, fasta."} +
                               basic_options_str +
                               basic_version_str;
        EXPECT_TRUE(ranges::equal((my_stdout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
    }
}

TEST(validator_test, output_file)
{
    test::tmp_filename tmp_name{"testbox.fasta"};
    test::tmp_filename tmp_name_2{"testbox_2.fasta"};
    test::tmp_filename tmp_name_3{"testbox_3.fa"};

    std::vector formats{std::string{"fa"}, std::string{"sam"}, std::string{"fasta"}};

    { // single file

        { // empty list of file.
            output_file_validator my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name.get_path()));
        }

        { // file does not exist.
            std::ofstream tmp_file_2(tmp_name_2.get_path());
            std::filesystem::path does_not_exists{tmp_name_2.get_path()};
            output_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exists), parser_invalid_argument);
        }

        { // file has wrong format.
            output_file_validator my_validator{std::vector{std::string{"sam"}}};
                EXPECT_THROW(my_validator(tmp_name.get_path()), parser_invalid_argument);
        }

        { // file has no extension.
            std::filesystem::path no_extension{tmp_name.get_path()};
            no_extension.replace_extension();
            output_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(no_extension), parser_invalid_argument);
        }

        std::filesystem::path file_out_path;

        // option
        std::filesystem::path path = tmp_name.get_path();
        const char * argv[] = {"./argument_parser_test", "-o", path.c_str()};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(file_out_path, 'o', "out-option", "desc",
                          option_spec::DEFAULT, output_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(file_out_path.string(), path.string());
    }

    { // file list.
        std::vector<std::filesystem::path> output_files;

        // option
        std::filesystem::path path = tmp_name.get_path();
        std::filesystem::path path_3 = tmp_name_3.get_path();

        const char * argv[] = {"./argument_parser_test", path.c_str(), path_3.c_str()};
        argument_parser parser("test_parser", 3, argv);
        parser.add_positional_option(output_files, "desc", output_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(output_files.size(), 2u);
        EXPECT_EQ(output_files[0].string(), path.string());
        EXPECT_EQ(output_files[1].string(), path_3.string());
    }

    // get help page message
    {
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        argument_parser parser("test_parser", 2, argv);
        parser.add_positional_option(path, "desc", output_file_validator{formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser"
                               "==========="
                               "POSITIONAL ARGUMENTS"
                               "    ARGUMENT-1 (std::filesystem::path)"
                               "          desc Default: \"\". Output file formats: fa, sam, fasta."} +
                               basic_options_str +
                               basic_version_str;
        EXPECT_TRUE(ranges::equal((my_stdout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
    }
}

TEST(validator_test, input_directory)
{
    test::tmp_filename tmp_name{"testbox.fasta"};

    { // directory

        { // has filename
            std::ofstream tmp_dir(tmp_name.get_path());
            input_directory_validator my_validator{};
            EXPECT_THROW(my_validator(tmp_name.get_path()), parser_invalid_argument);
        }

        { // read directory
            std::filesystem::path p = tmp_name.get_path();
            p.remove_filename();
            std::ofstream tmp_dir(p);
            input_directory_validator my_validator{};
            my_validator(p);
            EXPECT_NO_THROW(my_validator(p));

            std::filesystem::path dir_in_path;

            // option
            const char * argv[] = {"./argument_parser_test", "-i", p.c_str()};
            argument_parser parser("test_parser", 3, argv);
            parser.add_option(dir_in_path, 'i', "input-option", "desc",
                              option_spec::DEFAULT, input_directory_validator{});

            EXPECT_NO_THROW(parser.parse());
            EXPECT_EQ(p.string(), dir_in_path.string());
        }
    }

    {
        // get help page message
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        argument_parser parser("test_parser", 2, argv);
        parser.add_positional_option(path, "desc", input_directory_validator{});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser"
                               "==========="
                               "POSITIONAL ARGUMENTS"
                               "    ARGUMENT-1 (std::filesystem::path)"
                               "          desc Default: \"\". The input directory must exist and end with a valid directory separator."} +
                               basic_options_str +
                               basic_version_str;

        EXPECT_TRUE(ranges::equal((my_stdout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
    }
}

TEST(validator_test, output_directory)
{
    test::tmp_filename tmp_name{"testbox.fasta"};

    { // directory

        { // read directory
            std::filesystem::path p = tmp_name.get_path();
            p.remove_filename();
            output_directory_validator my_validator{};
            my_validator(p);
            EXPECT_NO_THROW();

            std::filesystem::path dir_out_path;

            // option
            const char * argv[] = {"./argument_parser_test", "-o", p.c_str()};
            argument_parser parser("test_parser", 3, argv);
            parser.add_option(dir_out_path, 'o', "output-option", "desc",
                              option_spec::DEFAULT, output_directory_validator{});

            EXPECT_NO_THROW(parser.parse());
            EXPECT_EQ(p.string(), dir_out_path.string());
        }

        { // has filename
            output_directory_validator my_validator{};
            EXPECT_THROW(my_validator(tmp_name.get_path()), parser_invalid_argument);
        }
    }

    {
        // get help page message
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        argument_parser parser("test_parser", 2, argv);
        parser.add_positional_option(path, "desc", output_directory_validator{});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser"
                               "==========="
                               "POSITIONAL ARGUMENTS"
                               "    ARGUMENT-1 (std::filesystem::path)"
                               "          desc Default: \"\". The output directory must not exist and end with a valid directory separator."} +
                               basic_options_str +
                               basic_version_str;

        EXPECT_TRUE(ranges::equal((my_stdout | std::view::filter(!is_space)), expected | std::view::filter(!is_space)));
    }
}

TEST(validator_test, arithmetic_range_validator_success)
{
    int option_value{0};
    std::vector<int> option_vector{};

    // option
    const char * argv[] = {"./argument_parser_test", "-i", "10"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 'i', "int-option", "desc",
                      option_spec::DEFAULT, arithmetic_range_validator(1, 20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // option - negative values
    const char * argv2[] = {"./argument_parser_test", "-i", "-10"};
    argument_parser parser2("test_parser", 3, argv2);
    parser2.add_option(option_value, 'i', "int-option", "desc",
                       option_spec::DEFAULT, arithmetic_range_validator(-20, 20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "10"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value, "desc", arithmetic_range_validator(1, 20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // positional option - negative values
    const char * argv4[] = {"./argument_parser_test", "--", "-10"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_value, "desc", arithmetic_range_validator(-20, 20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    argument_parser parser5("test_parser", 5, argv5);
    parser5.add_option(option_vector, 'i', "int-option", "desc",
                       option_spec::DEFAULT, arithmetic_range_validator(-50,50));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 48);

    // positional option - vector
    option_vector.clear();
    const char * argv6[] = {"./argument_parser_test", "--", "-10", "1"};
    argument_parser parser6("test_parser", 4, argv6);
    parser6.add_positional_option(option_vector, "desc", arithmetic_range_validator(-20,20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser6.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 1);

    // get help page message
    option_vector.clear();
    const char * argv7[] = {"./argument_parser_test", "-h"};
    argument_parser parser7("test_parser", 2, argv7);
    parser7.add_positional_option(option_vector, "desc", arithmetic_range_validator(-20,20));

    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser"
                           "==========="
                           "POSITIONAL ARGUMENTS"
                           "    ARGUMENT-1 (List of signed 32 bit integer's)"
                           "          desc Default: []. Value must be in range [-20,20]." +
                           basic_options_str +
                           basic_version_str);
    EXPECT_TRUE(ranges::equal((my_stdout   | std::view::filter(!is_space)),
                               expected | std::view::filter(!is_space)));

    // option - double value
    double double_option_value;
    const char * argv8[] = {"./argument_parser_test", "-i", "10.9"};
    argument_parser parser8("test_parser", 3, argv8);
    parser8.add_option(double_option_value, 'i', "double-option", "desc",
                       option_spec::DEFAULT, arithmetic_range_validator(1, 20));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser8.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_FLOAT_EQ(double_option_value, 10.9);
}

TEST(validator_test, arithmetic_range_validator_error)
{
    int option_value;
    std::vector<int> option_vector;

    // option - above max
    const char * argv[] = {"./argument_parser_test", "-i", "30"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, 'i', "int-option", "desc",
                      option_spec::DEFAULT, arithmetic_range_validator(1, 20));

    EXPECT_THROW(parser.parse(), validation_failed);

    // option - below min
    const char * argv2[] = {"./argument_parser_test", "-i", "-21"};
    argument_parser parser2("test_parser", 3, argv2);
    parser2.add_option(option_value, 'i', "int-option", "desc",
                       option_spec::DEFAULT, arithmetic_range_validator(-20, 20));

    EXPECT_THROW(parser2.parse(), validation_failed);

    // positional option - above max
    const char * argv3[] = {"./argument_parser_test", "30"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value, "desc", arithmetic_range_validator(1, 20));

    EXPECT_THROW(parser3.parse(), validation_failed);

    // positional option - below min
    const char * argv4[] = {"./argument_parser_test", "--", "-21"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_value, "desc", arithmetic_range_validator(-20, 20));

    EXPECT_THROW(parser4.parse(), validation_failed);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-100"};
    argument_parser parser5("test_parser", 3, argv5);
    parser5.add_option(option_vector, 'i', "int-option", "desc",
                       option_spec::DEFAULT, arithmetic_range_validator(-50, 50));

    EXPECT_THROW(parser5.parse(), validation_failed);

    // positional option - vector
    option_vector.clear();
    const char * argv6[] = {"./argument_parser_test", "--", "-10", "100"};
    argument_parser parser6("test_parser", 4, argv6);
    parser6.add_positional_option(option_vector, "desc", arithmetic_range_validator(-20, 20));

    EXPECT_THROW(parser6.parse(), validation_failed);

    // option - double value
    double double_option_value;
    const char * argv7[] = {"./argument_parser_test", "-i", "0.9"};
    argument_parser parser7("test_parser", 3, argv7);
    parser7.add_option(double_option_value, 'i', "double-option", "desc",
                       option_spec::DEFAULT, arithmetic_range_validator(1, 20));

    EXPECT_THROW(parser7.parse(), validation_failed);
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
                      option_spec::DEFAULT, value_list_validator({"ha", "ba", "ma"}));

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
    parser3.add_positional_option(option_value, "desc", value_list_validator({"ha", "ba", "ma"}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ma");

    // positional option - vector
    const char * argv4[] = {"./argument_parser_test", "ha", "ma"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_vector, "desc",
                                  value_list_validator({"ha", "ba", "ma"}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "ha");
    EXPECT_EQ(option_vector[1], "ma");

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    argument_parser parser5("test_parser", 5, argv5);
    parser5.add_option(option_vector_int, 'i', "int-option", "desc",
                       option_spec::DEFAULT, value_list_validator<int>({-10,48,50}));

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector_int[0], -10);
    EXPECT_EQ(option_vector_int[1], 48);

    // get help page message
    option_vector_int.clear();
    const char * argv7[] = {"./argument_parser_test", "-h"};
    argument_parser parser7("test_parser", 2, argv7);
    parser7.add_option(option_vector_int, 'i', "int-option", "desc",
                       option_spec::DEFAULT, value_list_validator<int>({-10,48,50}));

    option_vector_int.clear();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser"
                           "===========" +
                           basic_options_str +
                           "    -i, --int-option (List of signed 32 bit integer's)"
                           "          desc Default: []. Value must be one of [-10,48,50]." +
                           basic_version_str);
    EXPECT_TRUE(ranges::equal((my_stdout   | std::view::filter(!is_space)),
                               expected | std::view::filter(!is_space)));
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
                      option_spec::DEFAULT, value_list_validator({"ha", "ba", "ma"}));

    EXPECT_THROW(parser.parse(), validation_failed);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "30"};
    argument_parser parser3("test_parser", 2, argv3);
    parser3.add_positional_option(option_value_int, "desc", value_list_validator({0, 5, 10}));

    EXPECT_THROW(parser3.parse(), validation_failed);

    // positional option - vector
    const char * argv4[] = {"./argument_parser_test", "fo", "ma"};
    argument_parser parser4("test_parser", 3, argv4);
    parser4.add_positional_option(option_vector, "desc",
                                  value_list_validator({"ha", "ba", "ma"}));

    EXPECT_THROW(parser4.parse(), validation_failed);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "488"};
    argument_parser parser5("test_parser", 5, argv5);
    parser5.add_option(option_vector_int, 'i', "int-option", "desc",
                       option_spec::DEFAULT, value_list_validator<int>({-10,48,50}));

    EXPECT_THROW(parser5.parse(), validation_failed);
}

TEST(validator_test, regex_validator_success)
{
    std::string option_value;
    std::vector<std::string> option_vector;
    regex_validator email_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");
    regex_validator email_vector_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");

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
                                  regex_validator{"^chr[0-9]+"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "chr1");

    // positional option - vector
    const char * argv3[] = {"./argument_parser_test", "rollo", "bollo", "lollo"};
    argument_parser parser3("test_parser", 4, argv3);
    parser3.add_positional_option(option_vector, "desc",
                                  regex_validator{".*oll.*"});

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

    // get help page message
    option_vector.clear();
    const char * argv7[] = {"./argument_parser_test", "-h"};
    argument_parser parser7("test_parser", 2, argv7);
    parser7.add_option(option_vector, 's', "string-option", "desc",
                       option_spec::DEFAULT, email_vector_validator);

    option_vector.clear();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser"
                           "===========" +
                           basic_options_str +
                           "    -s, --string-option (List of std::string's)"
                           "          desc Default: []. Value must match the pattern '[a-zA-Z]+@[a-zA-Z]+\\.com'." +
                           basic_version_str);
    EXPECT_TRUE(ranges::equal((my_stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));
}

TEST(validator_test, regex_validator_error)
{
    std::string option_value;
    std::vector<std::string> option_vector;

    // option
    const char * argv[] = {"./argument_parser_test", "--string-option", "sally"};
    argument_parser parser("test_parser", 3, argv);
    parser.add_option(option_value, '\0', "string-option", "desc",
                      option_spec::DEFAULT, regex_validator("tt"));

    EXPECT_THROW(parser.parse(), validation_failed);

    // positional option
    const char * argv2[] = {"./argument_parser_test", "jessy"};
    argument_parser parser2("test_parser", 2, argv2);
    parser2.add_positional_option(option_value, "desc",
                                  regex_validator("[0-9]"));

    EXPECT_THROW(parser2.parse(), validation_failed);

    // positional option - vector
    const char * argv3[] = {"./argument_parser_test", "rollo", "bttllo", "lollo"};
    argument_parser parser3("test_parser", 4, argv3);
    parser3.add_positional_option(option_vector, "desc",
                                  regex_validator{".*oll.*"});

    EXPECT_THROW(parser3.parse(), validation_failed);

    // option - vector
    option_vector.clear();
    const char * argv4[] = {"./argument_parser_test", "-s", "gh", "-s", "tt"};
    argument_parser parser4("test_parser", 5, argv4);
    parser4.add_option(option_vector, 's', "", "desc",
                       option_spec::DEFAULT, regex_validator("tt"));

    EXPECT_THROW(parser4.parse(), validation_failed);
}

TEST(validator_test, chaining_validators)
{
    std::string option_value{};
    std::vector<std::string> option_vector{};
    regex_validator absolute_path_validator("(/[^/]+)+/.*\\.[^/\\.]+$");
    output_file_validator my_file_ext_validator({"sa", "so"});

    // option
    {
        const char * argv[] = {"./argument_parser_test", "-s", "/absolute/path/file.sa"};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(option_value, 's', "string-option", "desc",
                          option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, "/absolute/path/file.sa");
    }

    {
        const char * argv[] = {"./argument_parser_test", "-s", "relative/path/file.sa"};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(option_value, 's', "string-option", "desc",
                          option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);

        EXPECT_THROW(parser.parse(), validation_failed);
    }

    {
        const char * argv[] = {"./argument_parser_test", "-s", "/absoulte/path/file.notValidExtension"};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(option_value, 's', "string-option", "desc",
                          option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);

        EXPECT_THROW(parser.parse(), validation_failed);
    }

    // with temporary validators
    {
        const char * argv[] = {"./argument_parser_test", "-s", "/absolute/path/file.sa"};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(option_value, 's', "string-option", "desc",
                          option_spec::DEFAULT,
                          regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} |
                          output_file_validator{{"sa", "so"}});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, "/absolute/path/file.sa");
    }

    // three validators
    {
        const char * argv[] = {"./argument_parser_test", "-s", "/absolute/path/file.sa"};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(option_value, 's', "string-option", "desc",
                          option_spec::DEFAULT,
                          regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} |
                          output_file_validator{{"sa", "so"}} |
                          regex_validator{".*"});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, "/absolute/path/file.sa");
    }

    // help page message
    {
        const char * argv[] = {"./argument_parser_test", "-h"};
        argument_parser parser("test_parser", 2, argv);
        parser.add_option(option_value, 's', "string-option", "desc",
                          option_spec::DEFAULT,
                          regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} |
                          output_file_validator{{"sa", "so"}} |
                          regex_validator{".*"});

        testing::internal::CaptureStdout();
        option_value.clear();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser"
                               "===========" +
                               basic_options_str +
                               "    -s, --string-option (std::string)"
                               "          desc Default:. Value must match the pattern '(/[^/]+)+/.*\\.[^/\\.]+$'. "
                               "          Output file formats:  sa, so."
                               "          Value must match the pattern '.*'."} +
                               basic_version_str;
        EXPECT_TRUE(ranges::equal((my_stdout   | ranges::view::remove_if(is_space)),
                                   expected | ranges::view::remove_if(is_space)));
    }

    // chaining with a container option value type
    {
        std::vector<std::string> option_list_value{};
        const char * argv[] = {"./argument_parser_test", "-s", "/absolute/path/file.sa"};
        argument_parser parser("test_parser", 3, argv);
        parser.add_option(option_list_value, 's', "string-option", "desc",
                          option_spec::DEFAULT,
                          regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} | output_file_validator{{"sa", "so"}});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_list_value[0], "/absolute/path/file.sa");
    }
}
