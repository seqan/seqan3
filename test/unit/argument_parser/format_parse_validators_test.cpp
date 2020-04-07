// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/filesystem>
#include <seqan3/test/tmp_filename.hpp>

struct dummy_file
{

    struct format1
    {
        static inline std::vector<std::string> file_extensions{ {"fa"}, {"fasta"}};
    };

    struct format2
    {
        static inline std::vector<std::string> file_extensions{ {"sam"}, {"bam"}};
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
                                     "          Export the help page information. Value must be one of [html, man].\n"
                                     "\n"
                                     "  \n";

std::string const basic_version_str = "VERSION\n"
                                      "    Last update: \n"
                                      "    test_parser version: \n"
                                      "    SeqAn version: " + seqan3::seqan3_version + "\n";

namespace seqan3::detail
{
struct test_accessor
{
    static void set_terminal_width(seqan3::argument_parser & parser, unsigned terminal_width)
    {
        std::visit([terminal_width](auto & f)
        {
            if constexpr(std::is_same_v<decltype(f), seqan3::detail::format_help &>)
                f.layout = seqan3::detail::format_help::console_layout_struct{terminal_width};
        }, parser.format);
    }
};
} // seqan3::detail

using seqan3::detail::test_accessor;

TEST(validator_test, fullfill_concept)
{
    EXPECT_FALSE(seqan3::validator<int>);

    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<int>>);
    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<int> const>);
    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<int> &>);

    EXPECT_TRUE(seqan3::validator<seqan3::detail::default_validator<std::vector<int>>>);
    EXPECT_TRUE(seqan3::validator<seqan3::arithmetic_range_validator>);
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
    seqan3::test::tmp_filename tmp_name{"testbox.fasta"};
    seqan3::test::tmp_filename tmp_name_2{"testbox_2.fasta"};

    std::vector formats{std::string{"fa"}, std::string{"sam"}, std::string{"fasta"}};

    std::ofstream tmp_file(tmp_name.get_path());
    std::ofstream tmp_file_2(tmp_name_2.get_path());

    { // single file

        { // empty list of file.
            seqan3::input_file_validator my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name.get_path()));
        }

        { // file already exists.
            std::filesystem::path does_not_exist{tmp_name.get_path()};
            does_not_exist.replace_extension(".bam");
            seqan3::input_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exist), seqan3::validation_error);
        }

        { // file has wrong format.
            seqan3::input_file_validator my_validator{std::vector{std::string{"sam"}}};
            EXPECT_THROW(my_validator(tmp_name.get_path()), seqan3::validation_error);
        }

        { // file has no extension.
            std::filesystem::path does_not_exist{tmp_name.get_path()};
            does_not_exist.replace_extension();
            seqan3::input_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exist), seqan3::validation_error);
        }

        {  // read from file
            seqan3::input_file_validator<dummy_file> my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name.get_path()));
        }

        std::filesystem::path file_in_path;

        // option
        std::string const & path = tmp_name.get_path().string();
        const char * argv[] = {"./argument_parser_test", "-i", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(file_in_path, 'i', "int-option", "desc",
                          seqan3::option_spec::DEFAULT, seqan3::input_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(file_in_path.string(), path);
    }

    { // file list.
        std::vector<std::filesystem::path> input_files;

        // option
        std::string const & path = tmp_name.get_path().string();
        std::string const & path_2 = tmp_name_2.get_path().string();

        const char * argv[] = {"./argument_parser_test", path.c_str(), path_2.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(input_files, "desc", seqan3::input_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(input_files.size(), 2u);
        EXPECT_EQ(input_files[0].string(), path);
        EXPECT_EQ(input_files[1].string(), path_2);
    }

    { // get help page message
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(path, "desc", seqan3::input_file_validator{formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser\n"
                               "===========\n"
                               "\n"
                               "POSITIONAL ARGUMENTS\n"
                               "    ARGUMENT-1 (std::filesystem::path)\n"
                               "          desc The input file must exist and read permissions must be granted.\n"
                               "          Valid file extensions are: [fa, sam, fasta].\n"
                               "\n"} +
                               basic_options_str +
                               "\n" +
                               basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, input_file_ext_from_file)
{
    // Give as a template argument the seqan3 file type to get all valid extensions for this file.
    seqan3::input_file_validator<dummy_file> validator{};
    EXPECT_EQ(validator.get_help_page_message(), "The input file must exist and read permissions must be granted. "
                                                 "Valid file extensions are: [fa, fasta, sam, bam].");

    seqan3::input_file_validator validator2{};
    EXPECT_EQ(validator2.get_help_page_message(), "The input file must exist and read permissions must be granted.");
}

TEST(validator_test, output_file)
{
    seqan3::test::tmp_filename tmp_name{"testbox.fasta"};
    seqan3::test::tmp_filename tmp_name_2{"testbox_2.fasta"};
    seqan3::test::tmp_filename tmp_name_3{"testbox_3.fa"};

    std::vector formats{std::string{"fa"}, std::string{"sam"}, std::string{"fasta"}};

    { // single file

        { // empty list of file.
            seqan3::output_file_validator my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name.get_path()));
        }

        { // file does not exist.
            std::ofstream tmp_file_2(tmp_name_2.get_path());
            std::filesystem::path does_not_exist{tmp_name_2.get_path()};
            seqan3::output_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(does_not_exist), seqan3::validation_error);
        }

        { // file has wrong format.
            seqan3::output_file_validator my_validator{std::vector{std::string{"sam"}}};
            EXPECT_THROW(my_validator(tmp_name.get_path()), seqan3::validation_error);
        }

        { // file has no extension.
            std::filesystem::path no_extension{tmp_name.get_path()};
            no_extension.replace_extension();
            seqan3::output_file_validator my_validator{formats};
            EXPECT_THROW(my_validator(no_extension), seqan3::validation_error);
        }

        {  // read from file
            seqan3::output_file_validator<dummy_file> my_validator{};
            EXPECT_NO_THROW(my_validator(tmp_name.get_path()));
        }

        std::filesystem::path file_out_path;

        // option
        std::string const & path = tmp_name.get_path().string();
        const char * argv[] = {"./argument_parser_test", "-o", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(file_out_path, 'o', "out-option", "desc",
                          seqan3::option_spec::DEFAULT, seqan3::output_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(file_out_path.string(), path);
    }

    { // file list.
        std::vector<std::filesystem::path> output_files;

        // option
        std::string const & path = tmp_name.get_path().string();
        std::string const & path_3 = tmp_name_3.get_path().string();

        const char * argv[] = {"./argument_parser_test", path.c_str(), path_3.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(output_files, "desc", seqan3::output_file_validator{formats});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(output_files.size(), 2u);
        EXPECT_EQ(output_files[0].string(), path);
        EXPECT_EQ(output_files[1].string(), path_3);
    }

    // get help page message
    {
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_positional_option(path, "desc", seqan3::output_file_validator{formats});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser\n"
                               "===========\n"
                               "\n"
                               "POSITIONAL ARGUMENTS\n"
                               "    ARGUMENT-1 (std::filesystem::path)\n"
                               "          desc The output file must not exist already and write permissions\n"
                               "          must be granted. Valid file extensions are: [fa, sam, fasta].\n"
                               "\n"} +
                               basic_options_str +
                               "\n" +
                               basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, output_file_ext_from_file)
{
    // Give as a template argument the seqan3 file type to get all valid extensions for this file.
    seqan3::output_file_validator<dummy_file> validator{};
    EXPECT_EQ(validator.get_help_page_message(), "The output file must not exist already and write permissions must "
                                                 "be granted. Valid file extensions are: [fa, fasta, sam, bam].");

    seqan3::output_file_validator validator2{};
    EXPECT_EQ(validator2.get_help_page_message(), "The output file must not exist already and write permissions must "
                                                  "be granted.");
}

TEST(validator_test, input_directory)
{
    seqan3::test::tmp_filename tmp_name{"testbox.fasta"};

    { // directory

        { // has filename
            std::ofstream tmp_dir(tmp_name.get_path());
            seqan3::input_directory_validator my_validator{};
            EXPECT_THROW(my_validator(tmp_name.get_path()), seqan3::validation_error);
        }

        { // read directory
            std::filesystem::path p = tmp_name.get_path();
            p.remove_filename();
            std::ofstream tmp_dir(p);
            seqan3::input_directory_validator my_validator{};
            my_validator(p);
            EXPECT_NO_THROW(my_validator(p));

            std::filesystem::path dir_in_path;

            // option
            std::string const & path = p.string();
            const char * argv[] = {"./argument_parser_test", "-i", path.c_str()};
            seqan3::argument_parser parser{"test_parser", 3, argv, false};
            test_accessor::set_terminal_width(parser, 80);
            parser.add_option(dir_in_path, 'i', "input-option", "desc",
                              seqan3::option_spec::DEFAULT, seqan3::input_directory_validator{});

            EXPECT_NO_THROW(parser.parse());
            EXPECT_EQ(path, dir_in_path.string());
        }
    }

    {
        // get help page message
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, false};
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
                               "\n"} +
                               basic_options_str +
                               "\n" +
                               basic_version_str;

        EXPECT_EQ(my_stdout, expected);
    }
}

TEST(validator_test, output_directory)
{
    seqan3::test::tmp_filename tmp_name{"testbox.fasta"};

    { // read directory
        std::filesystem::path p = tmp_name.get_path();
        p.remove_filename();
        seqan3::output_directory_validator my_validator{};
        my_validator(p);
        EXPECT_NO_THROW();

        std::filesystem::path dir_out_path;

        // option
        std::string const & path = p.string();
        const char * argv[] = {"./argument_parser_test", "-o", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(dir_out_path, 'o', "output-option", "desc",
                          seqan3::option_spec::DEFAULT, seqan3::output_directory_validator{});

        EXPECT_NO_THROW(parser.parse());
        EXPECT_EQ(path, dir_out_path.string());
    }

    {
        // get help page message
        std::filesystem::path path;
        const char * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, false};
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
                                           "\n"} +
                                           basic_options_str +
                                           "\n" +
                                           basic_version_str;

        EXPECT_EQ(my_stdout, expected);
    }
}

#if __has_include(<filesystem>)
// Setting the permissions with perm_options is not available in the experimental/filesystem branch.

// In case this test is built as `root`, we want to exclude tests that check if certain missing permissions cause
// specific exceptions. For this, we check if read/write permissions are still available after the permissions were
// revoked. Note that `root` can always read/write even if user/group/all permissions are not set.

inline bool read_access(std::filesystem::path const & file)
{
    std::fstream stream;
    stream.open(file, std::ios::in);
    return !stream.fail();
}

inline bool write_access(std::filesystem::path const & file)
{
    if (std::filesystem::is_directory(file))
    {
        std::fstream stream;
        std::filesystem::path test_file{file};
        test_file /= "test";
        stream.open(test_file, std::ios::out);
        return !stream.fail();
    }
    else
    {
        std::fstream stream;
        stream.open(file, std::ios::out);
        return !stream.fail();
    }
}

TEST(validator_test, inputfile_not_readable)
{
    seqan3::test::tmp_filename tmp_name{"my_file.test"};
    std::filesystem::path tmp_file{tmp_name.get_path()};
    std::ofstream str{tmp_name.get_path()};

    EXPECT_NO_THROW(seqan3::input_file_validator{}(tmp_file));

    std::filesystem::permissions(tmp_file,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read |
                                 std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::remove);

    if (!read_access(tmp_file))
    {
        EXPECT_THROW(seqan3::input_file_validator{}(tmp_file), seqan3::validation_error);
    }

    std::filesystem::permissions(tmp_file,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read |
                                 std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::add);
}

TEST(validator_test, inputdir_not_readable)
{
    seqan3::test::tmp_filename tmp_name{"dir"};
    std::filesystem::path tmp_dir{tmp_name.get_path()};

    std::filesystem::create_directory(tmp_dir);

    EXPECT_NO_THROW(seqan3::input_directory_validator{}(tmp_dir));

    std::filesystem::permissions(tmp_dir,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read |
                                 std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::remove);

    if (!read_access(tmp_dir))
    {
        EXPECT_THROW(seqan3::input_directory_validator{}(tmp_dir), seqan3::validation_error);
    }

    std::filesystem::permissions(tmp_dir,
                                 std::filesystem::perms::owner_read | std::filesystem::perms::group_read |
                                 std::filesystem::perms::others_read,
                                 std::filesystem::perm_options::add);
}

TEST(validator_test, outputfile_not_writable)
{
    seqan3::test::tmp_filename tmp_name{"my_file.test"};
    std::filesystem::path tmp_file{tmp_name.get_path()};

    EXPECT_NO_THROW(seqan3::output_file_validator{}(tmp_file));

    // Parent path is not writable.
    std::filesystem::permissions(tmp_file.parent_path(),
                                 std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
                                 std::filesystem::perms::others_write,
                                 std::filesystem::perm_options::remove);

    if (!write_access(tmp_file))
    {
        EXPECT_THROW(seqan3::output_file_validator{}(tmp_file), seqan3::validation_error);
    }

    // make sure we can remove the directory.
    std::filesystem::permissions(tmp_file.parent_path(),
                                 std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
                                 std::filesystem::perms::others_write,
                                 std::filesystem::perm_options::add);
}

TEST(validator_test, outputdir_not_writable)
{
    { // parent dir is not writable.
        seqan3::test::tmp_filename tmp_name{"dir"};
        std::filesystem::path tmp_dir{tmp_name.get_path()};

        EXPECT_NO_THROW(seqan3::output_directory_validator{}(tmp_dir));
        EXPECT_FALSE(std::filesystem::exists(tmp_dir));
        // Parent path is not writable.
        std::filesystem::permissions(tmp_dir.parent_path(),
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
                                     std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::remove);

        if (!write_access(tmp_dir))
        {
            EXPECT_THROW(seqan3::output_directory_validator{}(tmp_dir), seqan3::validation_error);
        }

        // make sure we can remove the directory.
        std::filesystem::permissions(tmp_dir.parent_path(),
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
                                     std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::add);
    }

    {  // this dir is not writable
        seqan3::test::tmp_filename tmp_name{"dir"};
        std::filesystem::path tmp_dir{tmp_name.get_path()};

        std::filesystem::create_directory(tmp_dir);
        EXPECT_NO_THROW(seqan3::output_directory_validator{}(tmp_dir));

        // This path exists but is not writable.
        std::filesystem::permissions(tmp_dir,
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
                                     std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::remove);

        if (!write_access(tmp_dir))
        {
            EXPECT_THROW(seqan3::output_directory_validator{}(tmp_dir), seqan3::validation_error);
        }

        // make sure we can remove the directory.
        std::filesystem::permissions(tmp_dir,
                                     std::filesystem::perms::owner_write | std::filesystem::perms::group_write |
                                     std::filesystem::perms::others_write,
                                     std::filesystem::perm_options::add);
    }
}
#endif // __has_include(<filesystem>)

TEST(validator_test, arithmetic_range_validator_success)
{
    int option_value{0};
    std::vector<int> option_vector{};

    // option
    const char * argv[] = {"./argument_parser_test", "-i", "10"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value, 'i', "int-option", "desc",
                      seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // option - negative values
    const char * argv2[] = {"./argument_parser_test", "-i", "-10"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_option(option_value, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{-20, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "10"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{1, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, 10);

    // positional option - negative values
    const char * argv4[] = {"./argument_parser_test", "--", "-10"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, false};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{-20, 20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, -10);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, false};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{-50,50});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 48);

    // positional option - vector
    option_vector.clear();
    const char * argv6[] = {"./argument_parser_test", "--", "-10", "1"};
    seqan3::argument_parser parser6{"test_parser", 4, argv6, false};
    test_accessor::set_terminal_width(parser6, 80);
    parser6.add_positional_option(option_vector, "desc", seqan3::arithmetic_range_validator{-20,20});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser6.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], -10);
    EXPECT_EQ(option_vector[1], 1);

    // get help page message
    option_vector.clear();
    const char * argv7[] = {"./argument_parser_test", "-h"};
    seqan3::argument_parser parser7{"test_parser", 2, argv7, false};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_positional_option(option_vector, "desc", seqan3::arithmetic_range_validator{-20,20});

    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser\n"
                           "===========\n"
                           "\n"
                           "POSITIONAL ARGUMENTS\n"
                           "    ARGUMENT-1 (List of signed 32 bit integer's)\n"
                           "          desc Default: []. Value must be in range [-20,20].\n"
                           "\n" +
                           basic_options_str +
                           "\n" +
                           basic_version_str);
    EXPECT_EQ(my_stdout, expected);

    // option - double value
    double double_option_value;
    const char * argv8[] = {"./argument_parser_test", "-i", "10.9"};
    seqan3::argument_parser parser8{"test_parser", 3, argv8, false};
    test_accessor::set_terminal_width(parser8, 80);
    parser8.add_option(double_option_value, 'i', "double-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1, 20});

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
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value, 'i', "int-option", "desc",
                      seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1, 20});

    EXPECT_THROW(parser.parse(), seqan3::validation_error);

    // option - below min
    const char * argv2[] = {"./argument_parser_test", "-i", "-21"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_option(option_value, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{-20, 20});

    EXPECT_THROW(parser2.parse(), seqan3::validation_error);

    // positional option - above max
    const char * argv3[] = {"./argument_parser_test", "30"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{1, 20});

    EXPECT_THROW(parser3.parse(), seqan3::validation_error);

    // positional option - below min
    const char * argv4[] = {"./argument_parser_test", "--", "-21"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, false};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_value, "desc", seqan3::arithmetic_range_validator{-20, 20});

    EXPECT_THROW(parser4.parse(), seqan3::validation_error);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-100"};
    seqan3::argument_parser parser5{"test_parser", 3, argv5, false};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{-50, 50});

    EXPECT_THROW(parser5.parse(), seqan3::validation_error);

    // positional option - vector
    option_vector.clear();
    const char * argv6[] = {"./argument_parser_test", "--", "-10", "100"};
    seqan3::argument_parser parser6{"test_parser", 4, argv6, false};
    test_accessor::set_terminal_width(parser6, 80);
    parser6.add_positional_option(option_vector, "desc", seqan3::arithmetic_range_validator{-20, 20});

    EXPECT_THROW(parser6.parse(), seqan3::validation_error);

    // option - double value
    double double_option_value;
    const char * argv7[] = {"./argument_parser_test", "-i", "0.9"};
    seqan3::argument_parser parser7{"test_parser", 3, argv7, false};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_option(double_option_value, 'i', "double-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::arithmetic_range_validator{1, 20});

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
    // all arithmetic types are deduced to double in order to easily allow chaining of arithmetic validators
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<double>,
                 decltype(seqan3::value_list_validator{1})>));
    // except char
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<char>,
                 decltype(seqan3::value_list_validator{'c'})>));
    // The same holds for a range of arithmetic types
    std::vector v{1, 2, 3};
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<double>,
                 decltype(seqan3::value_list_validator{v})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<double>,
                 decltype(seqan3::value_list_validator{v | std::views::take(2)})>));
    std::vector v_char{'1', '2', '3'};
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<char>,
                 decltype(seqan3::value_list_validator{v_char})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<char>,
                 decltype(seqan3::value_list_validator{v_char | std::views::take(2)})>));
    // const char * is deduced to std::string
    std::vector v2{"ha", "ba", "ma"};
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>,
                 decltype(seqan3::value_list_validator{"ha"})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>,
                 decltype(seqan3::value_list_validator{"ha", "ba", "ma"})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>,
                 decltype(seqan3::value_list_validator{v2})>));
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<std::string>,
                 decltype(seqan3::value_list_validator{v2 | std::views::take(2)})>));
    // custom types are used as is
    EXPECT_TRUE((std::same_as<seqan3::value_list_validator<foo>,
                              decltype(seqan3::value_list_validator{foo::one, foo::two})>));

    // usage
    // -----
    std::string option_value;
    int option_value_int;
    std::vector<std::string> option_vector;
    std::vector<int> option_vector_int;

    // option
    std::vector<std::string> valid_str_values{"ha", "ba", "ma"};
    const char * argv[] = {"./argument_parser_test", "-s", "ba"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value, 's', "string-option", "desc",
                      seqan3::option_spec::DEFAULT,
                      seqan3::value_list_validator{valid_str_values | std::views::take(2)});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ba");

    // option with integers
    const char * argv2[] = {"./argument_parser_test", "-i", "-21"};
    seqan3::argument_parser parser2{"test_parser", 3, argv2, false};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_option(option_value_int, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::value_list_validator<int>{0, -21, 10});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value_int, -21);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "ma"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value, "desc", seqan3::value_list_validator{valid_str_values});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ma");

    // positional option - vector
    const char * argv4[] = {"./argument_parser_test", "ha", "ma"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, false};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_vector, "desc", seqan3::value_list_validator{"ha", "ba", "ma"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "ha");
    EXPECT_EQ(option_vector[1], "ma");

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "48"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, false};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector_int, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::value_list_validator<int>{-10, 48, 50});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser5.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector_int[0], -10);
    EXPECT_EQ(option_vector_int[1], 48);

    // get help page message
    option_vector_int.clear();
    const char * argv7[] = {"./argument_parser_test", "-h"};
    seqan3::argument_parser parser7{"test_parser", 2, argv7, false};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_option(option_vector_int, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::value_list_validator<int>{-10, 48, 50});

    option_vector_int.clear();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser\n"
                           "===========\n"
                           "\n" +
                           basic_options_str +
                           "    -i, --int-option (List of signed 32 bit integer's)\n"
                           "          desc Default: []. Value must be one of [-10,48,50].\n"
                           "\n" +
                           basic_version_str);
    EXPECT_EQ(my_stdout, expected);
}

TEST(validator_test, value_list_validator_error)
{
    std::string option_value;
    int option_value_int;
    std::vector<std::string> option_vector;
    std::vector<int> option_vector_int;

    // option
    const char * argv[] = {"./argument_parser_test", "-s", "sa"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value, 's', "string-option", "desc",
                      seqan3::option_spec::DEFAULT, seqan3::value_list_validator{"ha", "ba", "ma"});

    EXPECT_THROW(parser.parse(), seqan3::validation_error);

    // positional option
    const char * argv3[] = {"./argument_parser_test", "30"};
    seqan3::argument_parser parser3{"test_parser", 2, argv3, false};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_value_int, "desc", seqan3::value_list_validator{0, 5, 10});

    EXPECT_THROW(parser3.parse(), seqan3::validation_error);

    // positional option - vector
    const char * argv4[] = {"./argument_parser_test", "fo", "ma"};
    seqan3::argument_parser parser4{"test_parser", 3, argv4, false};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_positional_option(option_vector, "desc",
                                  seqan3::value_list_validator{"ha", "ba", "ma"});

    EXPECT_THROW(parser4.parse(), seqan3::validation_error);

    // option - vector
    const char * argv5[] = {"./argument_parser_test", "-i", "-10", "-i", "488"};
    seqan3::argument_parser parser5{"test_parser", 5, argv5, false};
    test_accessor::set_terminal_width(parser5, 80);
    parser5.add_option(option_vector_int, 'i', "int-option", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::value_list_validator<int>{-10, 48, 50});

    EXPECT_THROW(parser5.parse(), seqan3::validation_error);
}

TEST(validator_test, regex_validator_success)
{
    std::string option_value;
    std::vector<std::string> option_vector;
    seqan3::regex_validator email_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");
    seqan3::regex_validator email_vector_validator("[a-zA-Z]+@[a-zA-Z]+\\.com");

    // option
    const char * argv[] = {"./argument_parser_test", "-s", "ballo@rollo.com"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value, 's', "string-option", "desc",
                      seqan3::option_spec::DEFAULT, email_validator);

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "ballo@rollo.com");

    // positional option
    const char * argv2[] = {"./argument_parser_test", "chr1"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, false};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_positional_option(option_value, "desc",
                                  seqan3::regex_validator{"^chr[0-9]+"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser2.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_value, "chr1");

    // positional option - vector
    const char * argv3[] = {"./argument_parser_test", "rollo", "bollo", "lollo"};
    seqan3::argument_parser parser3{"test_parser", 4, argv3, false};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_vector, "desc",
                                  seqan3::regex_validator{".*oll.*"});

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser3.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "rollo");
    EXPECT_EQ(option_vector[1], "bollo");
    EXPECT_EQ(option_vector[2], "lollo");

    // option - vector
    option_vector.clear();
    const char * argv4[] = {"./argument_parser_test", "-s", "rita@rambo.com", "-s", "tina@rambo.com"};
    seqan3::argument_parser parser4{"test_parser", 5, argv4, false};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_option(option_vector, 's', "string-option", "desc",
                       seqan3::option_spec::DEFAULT, email_vector_validator);

    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser4.parse());
    EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
    EXPECT_EQ(option_vector[0], "rita@rambo.com");
    EXPECT_EQ(option_vector[1], "tina@rambo.com");

    // get help page message
    option_vector.clear();
    const char * argv7[] = {"./argument_parser_test", "-h"};
    seqan3::argument_parser parser7{"test_parser", 2, argv7, false};
    test_accessor::set_terminal_width(parser7, 80);
    parser7.add_option(option_vector, 's', "string-option", "desc",
                       seqan3::option_spec::DEFAULT, email_vector_validator);

    option_vector.clear();
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser7.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    std::string my_stdout = testing::internal::GetCapturedStdout();
    std::string expected = std::string("test_parser\n"
                           "===========\n"
                           "\n" +
                           basic_options_str +
                           "    -s, --string-option (List of std::string's)\n"
                           "          desc Default: []. Value must match the pattern\n"
                           "          '[a-zA-Z]+@[a-zA-Z]+\\.com'.\n"
                           "\n" +
                           basic_version_str);
    EXPECT_EQ(my_stdout, expected);
}

TEST(validator_test, regex_validator_error)
{
    std::string option_value;
    std::vector<std::string> option_vector;

    // option
    const char * argv[] = {"./argument_parser_test", "--string-option", "sally"};
    seqan3::argument_parser parser{"test_parser", 3, argv, false};
    test_accessor::set_terminal_width(parser, 80);
    parser.add_option(option_value, '\0', "string-option", "desc",
                      seqan3::option_spec::DEFAULT, seqan3::regex_validator{"tt"});

    EXPECT_THROW(parser.parse(), seqan3::validation_error);

    // positional option
    const char * argv2[] = {"./argument_parser_test", "jessy"};
    seqan3::argument_parser parser2{"test_parser", 2, argv2, false};
    test_accessor::set_terminal_width(parser2, 80);
    parser2.add_positional_option(option_value, "desc",
                                  seqan3::regex_validator{"[0-9]"});

    EXPECT_THROW(parser2.parse(), seqan3::validation_error);

    // positional option - vector
    const char * argv3[] = {"./argument_parser_test", "rollo", "bttllo", "lollo"};
    seqan3::argument_parser parser3{"test_parser", 4, argv3, false};
    test_accessor::set_terminal_width(parser3, 80);
    parser3.add_positional_option(option_vector, "desc",
                                  seqan3::regex_validator{".*oll.*"});

    EXPECT_THROW(parser3.parse(), seqan3::validation_error);

    // option - vector
    option_vector.clear();
    const char * argv4[] = {"./argument_parser_test", "-s", "gh", "-s", "tt"};
    seqan3::argument_parser parser4{"test_parser", 5, argv4, false};
    test_accessor::set_terminal_width(parser4, 80);
    parser4.add_option(option_vector, 's', "", "desc",
                       seqan3::option_spec::DEFAULT, seqan3::regex_validator{"tt"});

    EXPECT_THROW(parser4.parse(), seqan3::validation_error);
}

TEST(validator_test, chaining_validators)
{
    std::string option_value{};
    std::vector<std::string> option_vector{};
    seqan3::regex_validator absolute_path_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"};
    seqan3::output_file_validator my_file_ext_validator{{"sa", "so"}};

    seqan3::test::tmp_filename tmp_name{"file.sa"};
    std::filesystem::path invalid_extension{tmp_name.get_path()};
    invalid_extension.replace_extension(".invalid");

    // option
    {
        std::string const & path = tmp_name.get_path().string();
        const char * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, path);
    }

    {
        auto rel_path = tmp_name.get_path().relative_path().string();
        const char * argv[] = {"./argument_parser_test", "-s", rel_path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);

        EXPECT_THROW(parser.parse(), seqan3::validation_error);
    }

    {
        std::string const & path = invalid_extension.string();
        const char * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT, absolute_path_validator | my_file_ext_validator);

        EXPECT_THROW(parser.parse(), seqan3::validation_error);
    }

    // with temporary validators
    {
        std::string const & path = tmp_name.get_path().string();
        const char * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT,
                          seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} |
                          seqan3::output_file_validator{{"sa", "so"}});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, path);
    }

    // three validators
    {
        std::string const & path = tmp_name.get_path().string();
        const char * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT,
                          seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} |
                          seqan3::output_file_validator{{"sa", "so"}} |
                          seqan3::regex_validator{".*"});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_value, path);
    }

    // help page message
    {
        option_value.clear();
        const char * argv[] = {"./argument_parser_test", "-h"};
        seqan3::argument_parser parser{"test_parser", 2, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT,
                          seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} |
                          seqan3::output_file_validator{{"sa", "so"}} |
                          seqan3::regex_validator{".*"});

        testing::internal::CaptureStdout();
        EXPECT_EXIT(parser.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
        std::string my_stdout = testing::internal::GetCapturedStdout();
        std::string expected = std::string{"test_parser\n"
                               "===========\n"
                               "\n" +
                               basic_options_str +
                               "    -s, --string-option (std::string)\n"
                               "          desc Default: . Value must match the pattern '(/[^/]+)+/.*\\.[^/\\.]+$'.\n"
                               "          The output file must not exist already and write permissions must be\n"
                               "          granted. Valid file extensions are: [sa, so]. Value must match the\n"
                               "          pattern '.*'.\n"
                               "\n"} +
                               basic_version_str;
        EXPECT_EQ(my_stdout, expected);
    }

    // chaining with a container option value type
    {
        std::vector<std::string> option_list_value{};
        std::string const & path = tmp_name.get_path().string();
        const char * argv[] = {"./argument_parser_test", "-s", path.c_str()};
        seqan3::argument_parser parser{"test_parser", 3, argv, false};
        test_accessor::set_terminal_width(parser, 80);
        parser.add_option(option_list_value, 's', "string-option", "desc",
                          seqan3::option_spec::DEFAULT,
                          seqan3::regex_validator{"(/[^/]+)+/.*\\.[^/\\.]+$"} | seqan3::output_file_validator{{"sa",
                                                                                                               "so"}});

        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        EXPECT_TRUE((testing::internal::GetCapturedStderr()).empty());
        EXPECT_EQ(option_list_value[0], path);
    }
}
