// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <ranges>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/test/file_access.hpp>
#include <seqan3/test/tmp_directory.hpp>

struct dummy_file
{

    struct format1
    {
        static inline std::vector<std::string> file_extensions{{"fa"}, {"fasta"}};
    };

    struct format2
    {
        static inline std::vector<std::string> file_extensions{{"sam"}, {"bam"}};
    };

    using valid_formats = seqan3::type_list<format1, format2>;
};

std::string const basic_options_str = "OPTIONS\n"
                                      "\n"
                                      "  Basic options:\n"
                                      "    -h, --help\n"
                                      "          Prints the help page.\n"
                                      "    -hh, --advanced-help\n"
                                      "          Prints the help page including advanced options.\n"
                                      "    --version\n"
                                      "          Prints the version information.\n"
                                      "    --copyright\n"
                                      "          Prints the copyright/license information.\n"
                                      "    --export-help (std::string)\n"
                                      "          Export the help page information. Value must be one of [html, man].\n";

std::string const basic_version_str = "VERSION\n"
                                      "    Last update:\n"
                                      "    test_parser version:\n"
                                      "    SeqAn version: "
                                    + std::string{seqan3::seqan3_version_cstring} + "\n";

namespace seqan3::detail
{
struct test_accessor
{
    static void set_terminal_width(seqan3::argument_parser & parser, unsigned terminal_width)
    {
        std::visit(
            [terminal_width](auto & f)
            {
                if constexpr (std::is_same_v<decltype(f), seqan3::detail::format_help &>)
                    f.layout = seqan3::detail::format_help::console_layout_struct{terminal_width};
            },
            parser.format);
    }
};
} // namespace seqan3::detail

using seqan3::detail::test_accessor;

TEST(validator_test, fullfill_concept)
{
    EXPECT_FALSE(seqan3::validator<int>);

    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<int>>);
    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<int> const>);
    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<int> &>);

    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<std::vector<int>>>);
    EXPECT_TRUE(seqan3::validator<seqan3::arithmetic_range_validator<int>>);
    EXPECT_TRUE(seqan3::validator<seqan3::value_list_validator<double>>);
    EXPECT_TRUE(seqan3::validator<seqan3::value_list_validator<std::string>>);
    EXPECT_TRUE(seqan3::validator<seqan3::input_file_validator<>>);
    EXPECT_TRUE(seqan3::validator<seqan3::output_file_validator<>>);
    EXPECT_TRUE(seqan3::validator<seqan3::input_directory_validator>);
    EXPECT_TRUE(seqan3::validator<seqan3::output_directory_validator>);
    EXPECT_TRUE(seqan3::validator<seqan3::regex_validator>);

    EXPECT_TRUE(seqan3::validator<decltype(seqan3::input_file_validator{{"t"}} | seqan3::regex_validator{".*"})>);
}

TEST(validator_test, input_file)
{
    seqan3::test::tmp_directory tmp;
    auto tmp_name = tmp.path() / "testbox.fasta";
    auto tmp_name_2 = tmp.path() / "testbox_2.fasta";
    auto tmp_name_hidden = tmp.path() / ".testbox.fasta";
    auto tmp_name_multiple = tmp.path() / "testbox.fasta.txt";

    std::vector formats{std::string{"fa"}, std::string{"sam"}, std::string{"fasta"}, std::string{"fasta.txt"}};

    std::ofstream tmp_file(tmp_name);
    std::ofstream tmp_file_2(tmp_name_2);
    std::ofstream tmp_file_hidden(tmp_name_hidden);
    std::ofstream tmp_file_multiple(tmp_name_multiple);

    { // single file

        { // empty list of file.
            seqan3::input_file_validator my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name));
        }

        { // file already exists.
            std::filesystem::path does_not_exist{tmp_name};
            does_not_exist.replace_extension(".bam");
            seqan3::input_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exist), seqan3::validation_error);
        }

        { // file has wrong format.
            seqan3::input_file_validator my_validator{std::vector{std::string{"sam"}}};
            EXPECT_THROW(my_validator(tmp_name), seqan3::validation_error);
        }

        { // file has no extension.
            std::filesystem::path does_not_exist{tmp_name};
            does_not_exist.replace_extension();
            seqan3::input_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exist), seqan3::validation_error);
        }

        { // filename starts with dot.
            seqan3::input_file_validator my_validator{formats};
            EXPECT_NO_THROW(my_validator(tmp_name_hidden));
        }

        { // file has multiple extensions.
            seqan3::input_file_validator my_validator{formats};
            EXPECT_NO_THROW(my_validator(tmp_name_multiple));
        }

        { // read from file
            seqan3::input_file_validator<dummy_file> my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name));
        }

        std::filesystem::path file_in_path;

        // option
        std::string const & path = tmp_name.string();
        char const * argv[] = {"./argument_parser_test", "-i", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(file_in_path,
                          'i',
                          "int-option",
                          "desc",
                          seqan3::option_spec::standard,
                          seqan3::input_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(file_in_path.string(), path);
    }

    { // file list.
        std::vector<std::filesystem::path> input_files;

        // option
        std::string const & path = tmp_name.string();
        std::string const & path_2 = tmp_name_2.string();

        char const * argv[] = {"./argument_parser_test", path.c_str(), path_2.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(input_files, "desc", seqan3::input_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(input_files.size(), 2u);
        EXPECT_EQ(input_files[0].string(), path);
        EXPECT_EQ(input_files[1].string(), path_2);
    }

    { // get help page message
        std::filesystem::path path;
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(path, "desc", seqan3::input_file_validator{formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected =
            std::string{"test_parser\n"
                        "===========\n"
                        "\n"
                        "POSITIONAL ARGUMENTS\n"
                        "    ARGUMENT-1 (std::filesystem::path)\n"
                        "          desc The input file must exist and read permissions must be granted.\n"
                        "          Valid file extensions are: [fa, sam, fasta, fasta.txt].\n"
                        "\n"}
            + basic_options_str + "\n" + basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, input_file_ext_from_file)
{
    // Give as a template argument the seqan3 file type to get all valid extensions for this file.
    seqan3::input_file_validator<dummy_file> validator{};
    EXPECT_EQ(validator.get_help_page_message(),
              "The input file must exist and read permissions must be granted. "
              "Valid file extensions are: [fa, fasta, sam, bam].");

    seqan3::input_file_validator validator2{};
    EXPECT_EQ(validator2.get_help_page_message(), "The input file must exist and read permissions must be granted.");
}

TEST(validator_test, output_file)
{
    seqan3::test::tmp_directory tmp;
    auto tmp_name = tmp.path() / "testbox.fasta";
    auto not_existing_path = tmp_name;
    auto tmp_name_2 = tmp.path() / "testbox_2.fasta";
    std::ofstream tmp_file_2(tmp_name_2); // create file
    auto existing_path = tmp_name_2;
    auto tmp_name_3 = tmp.path() / "testbox_3.fa";
    auto hidden_name = tmp.path() / ".testbox.fasta";

    std::vector formats{std::string{"fa"}, std::string{"sam"}, std::string{"fasta"}, std::string{"fasta.txt"}};

    { // single file

        { // file does not exist (& no formats given)
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::open_or_create};
            EXPECT_NO_THROW(my_validator(not_existing_path));
            seqan3::output_file_validator my_validator2{seqan3::output_file_open_options::create_new};
            EXPECT_NO_THROW(my_validator2(not_existing_path));
            seqan3::output_file_validator my_validator3{}; // default: create_new
            EXPECT_NO_THROW(my_validator3(not_existing_path));
        }

        { // file does exist & overwriting is prohibited
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::create_new, formats};
            EXPECT_THROW(my_validator(existing_path), seqan3::validation_error);
        }

        { // file does exist but allow to overwrite it
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::open_or_create, formats};
            EXPECT_NO_THROW(my_validator(existing_path));
        }

        { // file has wrong format.
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::create_new,
                                                       std::vector{std::string{"sam"}}};
            EXPECT_THROW(my_validator(tmp_name), seqan3::validation_error);
        }

        { // file has no extension.
            std::filesystem::path no_extension{tmp_name};
            no_extension.replace_extension();
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::create_new, formats};
            EXPECT_THROW(my_validator(no_extension), seqan3::validation_error);
        }

        { // filename starts with dot.
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::create_new, formats};
            EXPECT_NO_THROW(my_validator(hidden_name));
        }

        { // file has multiple extensions.
            std::filesystem::path multiple_extension{tmp_name};
            multiple_extension.replace_extension("fasta.txt");
            seqan3::output_file_validator my_validator{seqan3::output_file_open_options::create_new, formats};
            EXPECT_NO_THROW(my_validator(multiple_extension));
        }

        { // read from file
            seqan3::output_file_validator<dummy_file> my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name));
        }

        std::filesystem::path file_out_path;

        // option
        std::string const & path = tmp_name.string();
        char const * argv[] = {"./argument_parser_test", "-o", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(file_out_path,
                          'o',
                          "out-option",
                          "desc",
                          seqan3::option_spec::standard,
                          seqan3::output_file_validator{seqan3::output_file_open_options::create_new, formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(file_out_path.string(), path);
    }

    { // file list.
        std::vector<std::filesystem::path> output_files;

        // option
        std::string const & path = tmp_name.string();
        std::string const & path_3 = tmp_name_3.string();

        char const * argv[] = {"./argument_parser_test", path.c_str(), path_3.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(
            output_files,
            "desc",
            seqan3::output_file_validator{seqan3::output_file_open_options::create_new, formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(output_files.size(), 2u);
        EXPECT_EQ(output_files[0].string(), path);
        EXPECT_EQ(output_files[1].string(), path_3);
    }

    // get help page message (overwriting prohibited)
    {
        std::filesystem::path path;
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(
            path,
            "desc",
            seqan3::output_file_validator{seqan3::output_file_open_options::create_new, formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected =
            std::string{"test_parser\n"
                        "===========\n"
                        "\n"
                        "POSITIONAL ARGUMENTS\n"
                        "    ARGUMENT-1 (std::filesystem::path)\n"
                        "          desc The output file must not exist already and write permissions\n"
                        "          must be granted. Valid file extensions are: [fa, sam, fasta,\n"
                        "          fasta.txt].\n"
                        "\n"}
            + basic_options_str + "\n" + basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }

    // get help page message (overwriting allowed)
    {
        std::filesystem::path path;
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(
            path,
            "desc",
            seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected =
            std::string{"test_parser\n"
                        "===========\n"
                        "\n"
                        "POSITIONAL ARGUMENTS\n"
                        "    ARGUMENT-1 (std::filesystem::path)\n"
                        "          desc Write permissions must be granted. Valid file extensions are:\n"
                        "          [fa, sam, fasta, fasta.txt].\n"
                        "\n"}
            + basic_options_str + "\n" + basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, output_file_ext_from_file)
{
    // Give as a template argument the seqan3 file type to get all valid extensions for this file.
    seqan3::output_file_validator<dummy_file> validator1{};
    EXPECT_EQ(validator1.get_help_page_message(),
              "The output file must not exist already and write permissions must "
              "be granted. Valid file extensions are: [fa, fasta, sam, bam].");

    seqan3::output_file_validator<dummy_file> validator2{seqan3::output_file_open_options::create_new};
    EXPECT_EQ(validator2.get_help_page_message(),
              "The output file must not exist already and write permissions must "
              "be granted. Valid file extensions are: [fa, fasta, sam, bam].");

    seqan3::output_file_validator<dummy_file> validator3{seqan3::output_file_open_options::open_or_create};
    EXPECT_EQ(validator3.get_help_page_message(),
              "Write permissions must be granted. Valid file extensions are: [fa, "
              "fasta, sam, bam].");

    seqan3::output_file_validator validator4{};
    EXPECT_EQ(validator4.get_help_page_message(),
              "The output file must not exist already and write permissions must "
              "be granted.");

    seqan3::output_file_validator validator5{seqan3::output_file_open_options::create_new};
    EXPECT_EQ(validator5.get_help_page_message(),
              "The output file must not exist already and write permissions must "
              "be granted.");

    seqan3::output_file_validator validator6{seqan3::output_file_open_options::open_or_create};
    EXPECT_EQ(validator6.get_help_page_message(), "Write permissions must be granted.");
}

TEST(validator_test, input_directory)
{
    seqan3::test::tmp_directory tmp{};
    auto tmp_name = tmp.path() / "testbox.fasta";

    { // directory

        { // has filename
            std::ofstream tmp_dir(tmp_name);
            seqan3::input_directory_validator my_validator{};
            EXPECT_THROW(my_validator(tmp_name), seqan3::validation_error);
        }

        { // read directory
            std::filesystem::path p = tmp_name;
            p.remove_filename();
            std::ofstream tmp_dir(p);
            seqan3::input_directory_validator my_validator{};
            my_validator(p);
            EXPECT_NO_THROW(my_validator(p));

            std::filesystem::path dir_in_path;

            // option
            std::string const & path = p.string();
            char const * argv[] = {"./argument_parser_test", "-i", path.c_str()};
            seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
            test_accessor::set_terminal_width(parser, 80);
            parser.add_option(dir_in_path,
                              'i',
                              "input-option",
                              "desc",
                              seqan3::option_spec::standard,
                              seqan3::input_directory_validator{});

            EXPECT_NO_THROW(parser.parse());
            EXPECT_EQ(path, dir_in_path.string());
        }
    }

    {
        // get help page message
        std::filesystem::path path;
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(path, "desc", seqan3::input_directory_validator{});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser\n"
                                           "===========\n"
                                           "\n"
                                           "POSITIONAL ARGUMENTS\n"
                                           "    ARGUMENT-1 (std::filesystem::path)\n"
                                           "          desc An existing, readable path for the input directory.\n"
                                           "\n"}
                             + basic_options_str + "\n" + basic_version_str;

        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, output_directory)
{
    seqan3::test::tmp_directory tmp{};

    { // read directory
        std::filesystem::path p = tmp.path() / "testbox.fasta";
        p.remove_filename();
        seqan3::output_directory_validator my_validator{};
        my_validator(p);
        EXPECT_NO_THROW();

        std::filesystem::path dir_out_path;

        // option
        std::string const & path = p.string();
        char const * argv[] = {"./argument_parser_test", "-o", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(dir_out_path,
                          'o',
                          "output-option",
                          "desc",
                          seqan3::option_spec::standard,
                          seqan3::output_directory_validator{});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(path, dir_out_path.string());
    }

    { // Parent path exists and is writable.
        std::filesystem::path tmp_child_name = tmp.path() / "dir/child_dir";
        std::filesystem::path tmp_child_dir{tmp_child_name};
        std::filesystem::path tmp_parent_path{tmp_child_dir.parent_path()};

        std::filesystem::create_directory(tmp_parent_path);

        EXPECT_TRUE(std::filesystem::exists(tmp_parent_path));
        EXPECT_NO_THROW(seqan3::output_directory_validator{}(tmp_child_dir));
    }

    {
        // get help page message
        std::filesystem::path path;
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(path, "desc", seqan3::output_directory_validator{});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser\n"
                                           "===========\n"
                                           "\n"
                                           "POSITIONAL ARGUMENTS\n"
                                           "    ARGUMENT-1 (std::filesystem::path)\n"
                                           "          desc A valid path for the output directory.\n"
                                           "\n"}
                             + basic_options_str + "\n" + basic_version_str;

        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, inputfile_not_readable)
{
    seqan3::test::tmp_directory tmp;
    auto tmp_name = tmp.path() / "my_file.test";
    std::filesystem::path tmp_file{tmp_name};
    std::ofstream str{tmp_name};

    EXPECT_NO_THROW(seqan3::input_file_validator{}(tmp_file));

    std::filesystem::permissions(tmp_file,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read
                                     | std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::remove);

    if (!seqan3::test::read_access(tmp_file)) // Do not execute with root permissions.
    {
        EXPECT_THROW(seqan3::input_file_validator{}(tmp_file), seqan3::validation_error);
    }

    std::filesystem::permissions(tmp_file,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read
                                     | std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::add);
}

TEST(validator_test, inputfile_not_regular)
{
    seqan3::test::tmp_directory tmp;
    std::filesystem::path filename = tmp.path() / "my_file.test";
    mkfifo(filename.c_str(), 0644);

    EXPECT_THROW(seqan3::input_file_validator{}(filename), seqan3::validation_error);
}

TEST(validator_test, inputdir_not_existing)
{
    seqan3::test::tmp_directory tmp;
    std::filesystem::path tmp_name = tmp.path() / "dir";

    std::filesystem::path not_existing_dir{tmp_name};

    EXPECT_THROW(seqan3::input_directory_validator{}(not_existing_dir), seqan3::validation_error);
}

TEST(validator_test, inputdir_not_readable)
{
    seqan3::test::tmp_directory tmp;
    auto tmp_dir = tmp.path() / "dir";

    std::filesystem::create_directory(tmp_dir);

    EXPECT_NO_THROW(seqan3::input_directory_validator{}(tmp_dir));

    std::filesystem::permissions(tmp_dir,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read
                                     | std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::remove);

    if (!seqan3::test::read_access(tmp_dir)) // Do not execute with root permissions.
    {
        EXPECT_THROW(seqan3::input_directory_validator{}(tmp_dir), seqan3::validation_error);
    }

    std::filesystem::permissions(tmp_dir,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read
                                     | std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::add);
}

TEST(validator_test, outputfile_not_writable)
{
    seqan3::test::tmp_directory tmp;
    auto tmp_file = tmp.path() / "my_file.test";

    EXPECT_NO_THROW(seqan3::output_file_validator{seqan3::output_file_open_options::create_new}(tmp_file));

    // Parent path is not writable.
    std::filesystem::permissions(tmp_file.parent_path(),
                                 std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                     | std::filesystem::perms::others_write,
                                 std::filesystem::perm_options::remove);

    if (!seqan3::test::write_access(tmp_file)) // Do not execute with root permissions.
    {
        EXPECT_THROW(seqan3::output_file_validator{seqan3::output_file_open_options::create_new}(tmp_file),
                     seqan3::validation_error);
    }

    // make sure we can remove the directory.
    std::filesystem::permissions(tmp_file.parent_path(),
                                 std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                     | std::filesystem::perms::others_write,
                                 std::filesystem::perm_options::add);
}

TEST(validator_test, outputdir_not_writable)
{
    { // parent dir is not writable.
        seqan3::test::tmp_directory tmp;
        auto tmp_dir = tmp.path() / "dir";

        EXPECT_NO_THROW(seqan3::output_file_validator{seqan3::output_file_open_options::create_new}(tmp_dir));
        EXPECT_FALSE(std::filesystem::exists(tmp_dir));

        // parent dir does not exist
        std::filesystem::path tmp_child_dir{tmp.path() / "dir/child_dir"};
        std::filesystem::path tmp_parent_dir{tmp_child_dir.parent_path()};

        EXPECT_THROW(seqan3::output_directory_validator{}(tmp_child_dir), seqan3::validation_error);

        // Directory exists but is not writable.
        std::filesystem::create_directory(tmp_dir);
        std::filesystem::permissions(tmp_dir,
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                         | std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::remove);

        EXPECT_TRUE(std::filesystem::exists(tmp_dir));
        if (!seqan3::test::write_access(tmp_dir)) // Do not execute with root permissions.
        {
            EXPECT_THROW(seqan3::output_directory_validator{}(tmp_dir), seqan3::validation_error);
        }

        // Parent path is not writable.
        std::filesystem::permissions(tmp_dir.parent_path(),
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                         | std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::remove);

        if (!seqan3::test::write_access(tmp_dir)) // Do not execute with root permissions.
        {
            EXPECT_THROW(seqan3::output_file_validator{seqan3::output_file_open_options::create_new}(tmp_dir),
                         seqan3::validation_error);
        }

        // make sure we can remove the directories.
        std::filesystem::permissions(tmp_dir,
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                         | std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::add);
        std::filesystem::permissions(tmp_dir.parent_path(),
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                         | std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::add);
    }

    { // this dir is not writable
        seqan3::test::tmp_directory tmp_dir{};

        std::filesystem::create_directory(tmp_dir.path());
        EXPECT_NO_THROW(seqan3::output_directory_validator{}(tmp_dir.path()));

        // This path exists but is not writable.
        std::filesystem::permissions(tmp_dir.path(),
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                         | std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::remove);

        if (!seqan3::test::write_access(tmp_dir.path())) // Do not execute with root permissions.
        {
            EXPECT_THROW(seqan3::output_file_validator{seqan3::output_file_open_options::create_new}(tmp_dir.path()),
                         seqan3::validation_error);
        }

        // make sure we can remove the directory.
        std::filesystem::permissions(tmp_dir.path(),
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write
                                         | std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::add);
    }
}

TEST(validator_test, arithmetic_range_validator_success)
{
    int option_value{0};
    std::vector<int> option_vector{};

    // option
    char const * argv[] = {"./argument_parser_test", "-i", "10"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value,
                      'i',
                      "int-option",
                      "desc",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // option - negative values
    char const * argv2[] = {"./argument_parser_test", "-i", "-10"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_option(option_value,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::arithmetic_range_validator{-20, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // positional option
    char const * argv3[] = {"./argument_parser_test", "10"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{1, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // positional option - negative values
    char const * argv4[] = {"./argument_parser_test", "--", "-10"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{-20, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // option - vector
    char const * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::arithmetic_range_validator{-50, 50});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 48);

    // positional option - vector
    option_vector.clear();
    char const * argv6[] = {"./argument_parser_test", "--", "-10", "1"};
    seqan3::argument_parser parser6{"test_parser", 4, argv6, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser6, 80);
    parser6.add_positional_option(option_vector, "desc", seqan3::arithmetic_range_validator{-20, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser6.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 1);

    // get help page message
    option_vector.clear();
    char const * argv7[] = {"./argument_parser_test", "-h"};
    seqan3::argument_parser parser7{"test_parser", 2, argv7, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_positional_option(option_vector, "desc", seqan3::arithmetic_range_validator{-20, 20});

    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser\n"
                                       "===========\n"
                                       "\n"
                                       "POSITIONAL ARGUMENTS\n"
                                       "    ARGUMENT-1 (List of signed 32 bit integer)\n"
                                       "          desc Default: []. Value must be in range [-20,20].\n"
                                       "\n"
                                       + basic_options_str + "\n" + basic_version_str);
    EXPECT_EQ(my_stdout, expected);

    // option - double value
    double double_option_value;
    char const * argv8[] = {"./argument_parser_test", "-i", "10.9"};
    seqan3::argument_parser parser8{"test_parser", 3, argv8, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser8, 80);
    parser8.add_option(double_option_value,
                       'i',
                       "double-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::arithmetic_range_validator{1, 20});

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
    char const * argv[] = {"./argument_parser_test", "-i", "30"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value,
                      'i',
                      "int-option",
                      "desc",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{1, 20});

    EXPECT_THROW(parser.parse(), seqan3::validation_error);

    // option - below min
    char const * argv2[] = {"./argument_parser_test", "-i", "-21"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_option(option_value,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::arithmetic_range_validator{-20, 20});

    EXPECT_THROW(parser2.parse(), seqan3::validation_error);

    // positional option - above max
    char const * argv3[] = {"./argument_parser_test", "30"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{1, 20});

    EXPECT_THROW(parser3.parse(), seqan3::validation_error);

    // positional option - below min
    char const * argv4[] = {"./argument_parser_test", "--", "-21"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{-20, 20});

    EXPECT_THROW(parser4.parse(), seqan3::validation_error);

    // option - vector
    char const * argv5[] = {"./argument_parser_test", "-i", "-100"};
    seqan3::argument_parser parser5{"test_parser", 3, argv5, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::arithmetic_range_validator{-50, 50});

    EXPECT_THROW(parser5.parse(), seqan3::validation_error);

    // positional option - vector
    option_vector.clear();
    char const * argv6[] = {"./argument_parser_test", "--", "-10", "100"};
    seqan3::argument_parser parser6{"test_parser", 4, argv6, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser6, 80);
    parser6.add_positional_option(option_vector, "desc", seqan3::arithmetic_range_validator{-20, 20});

    EXPECT_THROW(parser6.parse(), seqan3::validation_error);

    // option - double value
    double double_option_value;
    char const * argv7[] = {"./argument_parser_test", "-i", "0.9"};
    seqan3::argument_parser parser7{"test_parser", 3, argv7, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_option(double_option_value,
                       'i',
                       "double-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::arithmetic_range_validator{1, 20});

    EXPECT_THROW(parser7.parse(), seqan3::validation_error);
}

enum class foo
{
    one,
    two,
    three
};

TEST(validator_test, value_list_validator_success)
{
    // type deduction
    // --------------
    // all arithmetic types are deduced to their common type in order to easily allow chaining of arithmetic validators
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<int>, decltype(seqan3::value_list_validator{1})>));
    // except char
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<char>, decltype(seqan3::value_list_validator{'c'})>));
    // The same holds for a range of arithmetic types
    std::vector v{1, 2, 3};
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<int>, decltype(seqan3::value_list_validator{v})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<int>,
                              decltype(seqan3::value_list_validator{v | std::views::take(2)})>));
    std::vector v_char{'1', '2', '3'};
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<char>, decltype(seqan3::value_list_validator{v_char})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<char>,
                              decltype(seqan3::value_list_validator{v_char | std::views::take(2)})>));
    // const char * is deduced to std::string
    std::vector v2{"ha", "ba", "ma"};
    EXPECT_TRUE(
        (std::same_as<seqan3::value_list_validator<std::string>, decltype(seqan3::value_list_validator{"ha"})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>,
                              decltype(seqan3::value_list_validator{"ha", "ba", "ma"})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>, decltype(seqan3::value_list_validator{v2})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>,
                              decltype(seqan3::value_list_validator{v2 | std::views::take(2)})>));
    // custom types are used as is
    EXPECT_TRUE(
        (std::same_as<seqan3::value_list_validator<foo>, decltype(seqan3::value_list_validator{foo::one, foo::two})>));

    // usage
    // -----
    std::string option_value;
    int option_value_int;
    std::vector<std::string> option_vector;
    std::vector<int> option_vector_int;

    // option
    std::vector<std::string> valid_str_values{"ha", "ba", "ma"};
    char const * argv[] = {"./argument_parser_test", "-s", "ba"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value,
                      's',
                      "string-option",
                      "desc",
                      seqan3::option_spec::standard,
                      seqan3::value_list_validator{valid_str_values | std::views::take(2)});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ba");

    // option with integers
    char const * argv2[] = {"./argument_parser_test", "-i", "-21"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_option(option_value_int,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::value_list_validator<int>{0, -21, 10});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value_int, -21);

    // positional option
    char const * argv3[] = {"./argument_parser_test", "ma"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value, "desc", seqan3::value_list_validator{valid_str_values});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ma");

    // positional option - vector
    char const * argv4[] = {"./argument_parser_test", "ha", "ma"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_vector, "desc", seqan3::value_list_validator{"ha", "ba", "ma"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "ha");
    EXPECT_EQ(option_vector[1], "ma");

    // option - vector
    char const * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector_int,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::value_list_validator<int>{-10, 48, 50});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector_int[0], -10);
    EXPECT_EQ(option_vector_int[1], 48);

    // get help page message
    option_vector_int.clear();
    char const * argv7[] = {"./argument_parser_test", "-h"};
    seqan3::argument_parser parser7{"test_parser", 2, argv7, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_option(option_vector_int,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::value_list_validator<int>{-10, 48, 50});

    option_vector_int.clear();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser\n"
                                       "===========\n"
                                       "\n"
                                       + basic_options_str
                                       + "    -i, --int-option (List of signed 32 bit integer)\n"
                                         "          desc Default: []. Value must be one of [-10,48,50].\n\n"
                                       + basic_version_str);
    EXPECT_EQ(my_stdout, expected);
}

TEST(validator_test, value_list_validator_error)
{
    std::string option_value;
    int option_value_int;
    std::vector<std::string> option_vector;
    std::vector<int> option_vector_int;

    // option
    char const * argv[] = {"./argument_parser_test", "-s", "sa"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value,
                      's',
                      "string-option",
                      "desc",
                      seqan3::option_spec::standard,
                      seqan3::value_list_validator{"ha", "ba", "ma"});

    EXPECT_THROW(parser.parse(), seqan3::validation_error);

    // positional option
    char const * argv3[] = {"./argument_parser_test", "30"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value_int, "desc", seqan3::value_list_validator{0, 5, 10});

    EXPECT_THROW(parser3.parse(), seqan3::validation_error);

    // positional option - vector
    char const * argv4[] = {"./argument_parser_test", "fo", "ma"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_vector, "desc", seqan3::value_list_validator{"ha", "ba", "ma"});

    EXPECT_THROW(parser4.parse(), seqan3::validation_error);

    // option - vector
    char const * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "488"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector_int,
                       'i',
                       "int-option",
                       "desc",
                       seqan3::option_spec::standard,
                       seqan3::value_list_validator<int>{-10, 48, 50});

    EXPECT_THROW(parser5.parse(), seqan3::validation_error);
}

TEST(validator_test, regex_validator_success)
{
    std::string option_value;
    std::vector<std::string> option_vector;
    seqan3::regex_validator email_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");
    seqan3::regex_validator email_vector_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");

    { // option
        char const * argv[] = {"./argument_parser_test", "-s", "ballo@rollo.com"};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc", seqan3::option_spec::standard, email_validator);

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, "ballo@rollo.com");
    }

    { // positional option
        char const * argv[] = {"./argument_parser_test", "chr1"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(option_value, "desc", seqan3::regex_validator{"^chr[0-9]+"});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, "chr1");
    }

    { // positional option - vector
        char const * argv[] = {"./argument_parser_test", "rollo", "bollo", "lollo"};
        seqan3::argument_parser parser{"test_parser", 4, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(option_vector, "desc", seqan3::regex_validator{".*oll.*"});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_vector[0], "rollo");
        EXPECT_EQ(option_vector[1], "bollo");
        EXPECT_EQ(option_vector[2], "lollo");
    }

    { // option - vector
        option_vector.clear();
        char const * argv[] = {"./argument_parser_test", "-s", "rita@rambo.com", "-s", "tina@rambo.com"};
        seqan3::argument_parser parser{"test_parser", 5, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_vector,
                          's',
                          "string-option",
                          "desc",
                          seqan3::option_spec::standard,
                          email_vector_validator);

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_vector[0], "rita@rambo.com");
        EXPECT_EQ(option_vector[1], "tina@rambo.com");
    }

    { // option - std::filesystem::path
        std::filesystem::path path_option;
        char const * argv[] = {"./argument_parser_test", "-s", "rita@rambo.com"};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(path_option,
                          's',
                          "string-option",
                          "desc",
                          seqan3::option_spec::standard,
                          email_vector_validator);

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(path_option, "rita@rambo.com");
    }

    { // get help page message
        option_vector.clear();
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_vector,
                          's',
                          "string-option",
                          "desc",
                          seqan3::option_spec::standard,
                          email_vector_validator);

        option_vector.clear();
        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string("test_parser\n"
                                           "===========\n"
                                           "\n"
                                           + basic_options_str
                                           + "    -s, --string-option (List of std::string)\n"
                                             "          desc Default: []. Value must match the pattern\n"
                                             "          '[a-zA-Z]+@[a-zA-Z]+\\.com'.\n"
                                             "\n"
                                           + basic_version_str);
        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, regex_validator_error)
{
    std::string option_value;
    std::vector<std::string> option_vector;

    // option
    char const * argv[] = {"./argument_parser_test", "--string-option", "sally"};
    seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value,
                      '\0',
                      "string-option",
                      "desc",
                      seqan3::option_spec::standard,
                      seqan3::regex_validator{"tt"});

    EXPECT_THROW(parser.parse(), seqan3::validation_error);

    // positional option
    char const * argv2[] = {"./argument_parser_test", "jessy"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_positional_option(option_value, "desc", seqan3::regex_validator{"[0-9]"});

    EXPECT_THROW(parser2.parse(), seqan3::validation_error);

    // positional option - vector
    char const * argv3[] = {"./argument_parser_test", "rollo", "bttllo", "lollo"};
    seqan3::argument_parser parser3{"test_parser", 4, argv3, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_vector, "desc", seqan3::regex_validator{".*oll.*"});

    EXPECT_THROW(parser3.parse(), seqan3::validation_error);

    // option - vector
    option_vector.clear();
    char const * argv4[] = {"./argument_parser_test", "-s", "gh", "-s", "tt"};
    seqan3::argument_parser parser4{"test_parser", 5, argv4, seqan3::update_notifications::off};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_option(option_vector, 's', "", "desc", seqan3::option_spec::standard, seqan3::regex_validator{"tt"});

    EXPECT_THROW(parser4.parse(), seqan3::validation_error);
}

TEST(validator_test, chaining_validators_common_type)
{
    // chaining integral options stay integral
    {
        int max_int = std::numeric_limits<int>::max();
        std::vector v_int{1, 2, 3, max_int};
        std::vector v_unsigned{4u, static_cast<unsigned>(max_int)};

        EXPECT_TRUE((std::same_as<std::vector<int>, decltype(v_int)>));
        EXPECT_TRUE((std::same_as<std::vector<unsigned>, decltype(v_unsigned)>));

        seqan3::value_list_validator validator_int{v_int};
        seqan3::value_list_validator validator_unsigned{v_unsigned};

        EXPECT_TRUE((std::same_as<seqan3::value_list_validator<int>, decltype(validator_int)>));
        EXPECT_TRUE((std::same_as<int, decltype(validator_int)::option_value_type>));

        EXPECT_TRUE((std::same_as<seqan3::value_list_validator<unsigned>, decltype(validator_unsigned)>));
        EXPECT_TRUE((std::same_as<unsigned, decltype(validator_unsigned)::option_value_type>));

        auto validator = validator_int | validator_unsigned;

        EXPECT_TRUE((std::same_as<unsigned, std::common_type_t<int, unsigned>>));
        EXPECT_TRUE((std::same_as<unsigned, decltype(validator)::option_value_type>));

        // max_int is part of both validators
        EXPECT_NO_THROW(validator_int(max_int));
        EXPECT_NO_THROW(validator_unsigned(max_int));
        EXPECT_NO_THROW(validator(max_int));
    }

    // chaining mixed arithmetic options will be highest common arithmetic type
    {
        // note: this number is not representable by double and multiple integer values represent the same double value
        int64_t max_int64 = std::numeric_limits<int64_t>::max();
        EXPECT_EQ(static_cast<double>(max_int64), static_cast<double>(max_int64 - 1));

        std::vector<int64_t> v_int64{1, 2, 3, max_int64};
        std::vector<uint64_t> v_uint64{4u, static_cast<uint64_t>(max_int64)};
        std::vector<double> v_double{4.0, static_cast<double>(max_int64)};

        seqan3::value_list_validator validator_int64{v_int64};
        seqan3::value_list_validator validator_uint64{v_uint64};
        seqan3::value_list_validator validator_double{v_double};

        EXPECT_TRUE((std::same_as<seqan3::value_list_validator<int64_t>, decltype(validator_int64)>));
        EXPECT_TRUE((std::same_as<int64_t, decltype(validator_int64)::option_value_type>));

        EXPECT_TRUE((std::same_as<seqan3::value_list_validator<uint64_t>, decltype(validator_uint64)>));
        EXPECT_TRUE((std::same_as<uint64_t, decltype(validator_uint64)::option_value_type>));

        EXPECT_TRUE((std::same_as<seqan3::value_list_validator<double>, decltype(validator_double)>));
        EXPECT_TRUE((std::same_as<double, decltype(validator_double)::option_value_type>));

        auto validator = validator_int64 | validator_uint64 | validator_double;

        EXPECT_TRUE((std::same_as<double, std::common_type_t<int64_t, uint64_t, double>>));
        EXPECT_TRUE((std::same_as<double, decltype(validator)::option_value_type>));

        // max_int64 is an exact match for the two integral validators
        // note: double will decay the integer to a smaller double value,
        //       but this is consistent, because it is the same given value
        // note: chained validator passes the value as it is through,
        //       so validator_[u]int64 will be called with the integer value
        EXPECT_NO_THROW(validator_int64(max_int64));
        EXPECT_NO_THROW(validator_uint64(max_int64));
        EXPECT_NO_THROW(validator_double(max_int64));
        EXPECT_NO_THROW(validator(max_int64));

        // integers have exact match
        // note: double accepts that value, even though it is not within that list.
        EXPECT_THROW(validator_int64(max_int64 - 1), seqan3::validation_error);
        EXPECT_THROW(validator_uint64(max_int64 - 1), seqan3::validation_error);
        EXPECT_NO_THROW(validator_double(max_int64 - 1));
        EXPECT_THROW(validator(max_int64 - 1), seqan3::validation_error);
    }
}

TEST(validator_test, chaining_validators)
{
    std::string option_value{};
    std::vector<std::string> option_vector{};
    seqan3::regex_validator absolute_path_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"};
    seqan3::output_file_validator my_file_ext_validator{seqan3::output_file_open_options::create_new, {"sa", "so"}};

    seqan3::test::tmp_directory tmp;
    auto tmp_name = tmp.path() / "file.sa";

    std::filesystem::path invalid_extension{tmp_name};
    invalid_extension.replace_extension(".invalid");

    // option
    {
        std::string const & path = tmp_name.string();
        char const * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value,
                          's',
                          "string-option",
                          "desc",
                          seqan3::option_spec::standard,
                          absolute_path_validator | my_file_ext_validator);

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, path);
    }

    {
        auto rel_path = tmp_name.relative_path().string();
        char const * argv[] = {"./argument_parser_test", "-s", rel_path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value,
                          's',
                          "string-option",
                          "desc",
                          seqan3::option_spec::standard,
                          absolute_path_validator | my_file_ext_validator);

        EXPECT_THROW(parser.parse(), seqan3::validation_error);
    }

    {
        std::string const & path = invalid_extension.string();
        char const * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value,
                          's',
                          "string-option",
                          "desc",
                          seqan3::option_spec::standard,
                          absolute_path_validator | my_file_ext_validator);

        EXPECT_THROW(parser.parse(), seqan3::validation_error);
    }

    // with temporary validators
    {
        std::string const & path = tmp_name.string();
        char const * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(
            option_value,
            's',
            "string-option",
            "desc",
            seqan3::option_spec::standard,
            seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"}
                | seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sa", "so"}});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, path);
    }

    // three validators
    {
        std::string const & path = tmp_name.string();
        char const * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(
            option_value,
            's',
            "string-option",
            "desc",
            seqan3::option_spec::standard,
            seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"}
                | seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sa", "so"}}
                | seqan3::regex_validator{".*"});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, path);
    }

    // help page message
    {
        option_value.clear();
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(
            option_value,
            's',
            "string-option",
            "desc",
            seqan3::option_spec::standard,
            seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"}
                | seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sa", "so"}}
                | seqan3::regex_validator{".*"});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected =
            std::string{"test_parser\n"
                        "===========\n"
                        "\n"
                        + basic_options_str
                        + "    -s, --string-option (std::string)\n"
                          "          desc Default: . Value must match the pattern '(/[^/]+)+/.*\\.[^/\\.]+$'.\n"
                          "          The output file must not exist already and write permissions must be\n"
                          "          granted. Valid file extensions are: [sa, so]. Value must match the\n"
                          "          pattern '.*'.\n"
                          "\n"}
            + basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }

    // help page message (allow overwriting)
    {
        option_value.clear();
        char const * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(
            option_value,
            's',
            "string-option",
            "desc",
            seqan3::option_spec::standard,
            seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"}
                | seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"sa", "so"}}
                | seqan3::regex_validator{".*"});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected =
            std::string{"test_parser\n"
                        "===========\n"
                        "\n"
                        + basic_options_str
                        + "    -s, --string-option (std::string)\n"
                          "          desc Default: . Value must match the pattern '(/[^/]+)+/.*\\.[^/\\.]+$'.\n"
                          "          Write permissions must be granted. Valid file extensions are: [sa,\n"
                          "          so]. Value must match the pattern '.*'.\n"
                          "\n"}
            + basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }

    // chaining with a container option value type
    {
        std::vector<std::string> option_list_value{};
        std::string const & path = tmp_name.string();
        char const * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, seqan3::update_notifications::off};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(
            option_list_value,
            's',
            "string-option",
            "desc",
            seqan3::option_spec::standard,
            seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"}
                | seqan3::output_file_validator{seqan3::output_file_open_options::create_new, {"sa", "so"}});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_list_value[0], path);
    }
}
