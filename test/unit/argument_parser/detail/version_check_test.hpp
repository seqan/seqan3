// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <stdlib.h>

#include <seqan3/argument_parser/all.hpp>

using namespace seqan3;

//------------------------------------------------------------------------------
// test fixtures
//------------------------------------------------------------------------------

struct version_check : public ::testing::Test
{
    char const * const OPTION_VERSION_CHECK = "--version-check";
    char const * const OPTION_OFF = "0";
    char const * const OPTION_ON = "1";

    std::chrono::duration<long int>::rep const TIME_NOW =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    std::string const APP_NAME = "test_version_check_" + app_name_component + "_" + std::to_string(TIME_NOW); // avoid name conflicts.

    std::filesystem::path const PATH = detail::version_checker::get_path();

    std::filesystem::path const APP_VERSION_FILENAME = PATH / (APP_NAME + ".version");

#if defined(NDEBUG)
    std::filesystem::path const APP_TIMESTAMP_FILENAME = PATH / (APP_NAME + "_usr.timestamp");
#else
    std::filesystem::path const APP_TIMESTAMP_FILENAME = PATH / (APP_NAME + "_dev.timestamp");
#endif // defined(NDEBUG)

    std::string const SEQAN_VERSION = std::to_string(SEQAN3_VERSION_MAJOR) + "." +
                                      std::to_string(SEQAN3_VERSION_MINOR) + "." +
                                      std::to_string(SEQAN3_VERSION_PATCH);

    std::regex timestamp_regex{"^[[:digit:]]+$"}; // only digits

    std::tuple<std::string, std::string, bool> simulate_argument_parser(int argc, const char ** argv)
    {
        // make sure that the environment variable is not set
        char * env{std::getenv("SEQAN3_NO_VERSION_CHECK")};
        if (env != nullptr)
            unsetenv("SEQAN3_NO_VERSION_CHECK");

        bool app_call_succeeded{false};

        argument_parser parser{APP_NAME, argc, argv};
        parser.info.version = "2.3.4";

        // In case we don't want to specify --version-check but avoid that short help format will be set (no arguments)
        bool dummy{};
        parser.add_flag(dummy, 'f', "dummy-flag", "A dummy flag.");

        testing::internal::CaptureStdout();
        testing::internal::CaptureStderr();
        EXPECT_NO_THROW(parser.parse());
        std::string out = testing::internal::GetCapturedStdout();
        std::string err = testing::internal::GetCapturedStderr();

        // call future.get() to artificially wait for the thread to finish and avoid
        // any interference with following tests
        if (parser.version_check_future.valid())
            app_call_succeeded = parser.version_check_future.get();

        if (env != nullptr)
            setenv("SEQAN3_NO_VERSION_CHECK", env, 1);

        return {out, err, app_call_succeeded};
    }

    bool remove_files_from_path()
    {
        return (!std::filesystem::exists(APP_VERSION_FILENAME)   ||
                 std::filesystem::remove(APP_VERSION_FILENAME))  &&
               (!std::filesystem::exists(APP_TIMESTAMP_FILENAME) ||
                 std::filesystem::remove(APP_TIMESTAMP_FILENAME));
    }

    template <typename message_type>
    bool create_file(std::filesystem::path const & filename, message_type const & message)
    {
        std::ofstream out_file{filename};

        if (!out_file.is_open())
            return false;

        out_file << message;
        out_file.close();

        return true;
    }

    std::string read_file(std::filesystem::path const & filename)
    {
        std::ifstream in_file{filename};
        std::string line{};

        if (in_file.is_open())
        {
            std::getline(in_file, line);
            in_file.close();
        }

        return line;
    }
};

struct sanity_checks : public version_check
{};

//------------------------------------------------------------------------------
// sanity checks
//------------------------------------------------------------------------------

// even if the homedir might not be writable at least the tmp dir should be
TEST_F(sanity_checks, path_availability)
{
    EXPECT_FALSE(PATH.empty()) << "No writable directory found. All other tests cannot be trusted!";
}

TEST_F(sanity_checks, create_and_delete_files)
{
    EXPECT_TRUE(create_file(APP_VERSION_FILENAME, "20.5.9"));
    EXPECT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW));

    EXPECT_TRUE(std::filesystem::exists(APP_VERSION_FILENAME));
    EXPECT_TRUE(std::filesystem::exists(APP_TIMESTAMP_FILENAME));

    EXPECT_TRUE(remove_files_from_path()); // clear files again
    EXPECT_FALSE(std::filesystem::exists(APP_VERSION_FILENAME));
    EXPECT_FALSE(std::filesystem::exists(APP_TIMESTAMP_FILENAME));
}

//------------------------------------------------------------------------------
// version checks
//------------------------------------------------------------------------------

TEST_F(version_check, option_on)
{
    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_ON};

    auto [out, err, app_call_succeeded] = simulate_argument_parser(3, argv);

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    // no timestamp is written since the decision was made explicitly
    if (app_call_succeeded)
    {
        EXPECT_TRUE(std::filesystem::exists(APP_VERSION_FILENAME));
    }
    else
    {
        std::cout << "App call did not succeed (server offline?) and could thus not be tested." << std::endl;
    }

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

// Note that we cannot test interactiveness because google test captures std::cin and thus
// seqan3::detail::is_terminal() is always false
TEST_F(version_check, option_implicitely_on)
{
    const char * argv[2] = {APP_NAME.c_str(), "-f"};

    auto [out, err, app_call_succeeded] = simulate_argument_parser(2, argv);

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "\n#######################################################################\n"
                   "   Automatic Update Notifications\n"
                   "#######################################################################\n"
                   " This app performs automatic checks for updates. For more information\n"
                   " see: https://github.com/seqan/seqan3/wiki/Update-Notifications\n"
                   "#######################################################################\n\n");

    // make sure that all files now exist
    EXPECT_TRUE(std::filesystem::exists(APP_TIMESTAMP_FILENAME)) << APP_TIMESTAMP_FILENAME;
    EXPECT_TRUE(std::regex_match(read_file(APP_TIMESTAMP_FILENAME), timestamp_regex));

    if (app_call_succeeded)
    {
        EXPECT_TRUE(std::filesystem::exists(APP_VERSION_FILENAME));
    }
    else
    {
        std::cout << "App call did not succeed (server offline?) and could thus not be tested." << std::endl;
    }

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, time_out) // while implicitly on
{
     const char * argv[2] = {APP_NAME.c_str(), "-f"};

    // create timestamp files
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW));

    auto [out, err, app_call_succeeded] = simulate_argument_parser(2, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    EXPECT_FALSE(std::filesystem::exists(APP_VERSION_FILENAME));

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, environment_variable_set)
{
    // store variable for resetting it
    char * env{std::getenv("SEQAN3_NO_VERSION_CHECK")};
    setenv("SEQAN3_NO_VERSION_CHECK", "foo", 1);

    const char * argv[2] = {APP_NAME.c_str(), "-f"};

    argument_parser parser{APP_NAME, 2, argv};
    parser.info.version = "2.3.4";
    bool dummy{};
    parser.add_flag(dummy, 'f', "dummy-flag", "A dummy flag.");

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    std::string out = testing::internal::GetCapturedStdout();
    std::string err = testing::internal::GetCapturedStderr();

    if (parser.version_check_future.valid())
        parser.version_check_future.get();

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    // if environment variable is set, no cookies are written
    EXPECT_FALSE(std::filesystem::exists(APP_TIMESTAMP_FILENAME)) << APP_TIMESTAMP_FILENAME;
    EXPECT_FALSE(std::filesystem::exists(APP_VERSION_FILENAME));

    if (env == nullptr)
        unsetenv("SEQAN3_NO_VERSION_CHECK");
    else
        setenv("SEQAN3_NO_VERSION_CHECK", env, 1);

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, option_off)
{
    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_OFF};

    auto [out, err, app_call_succeeded] = simulate_argument_parser(3, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    // no timestamp is written since the decision was made explicitly
    EXPECT_FALSE(std::filesystem::exists(APP_VERSION_FILENAME)) << APP_VERSION_FILENAME;

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

// case: the current argument parser has a smaller seqan version than is present in the version file
#if !defined(NDEBUG)
TEST_F(version_check, smaller_SEQAN_VERSION)
{
    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_ON};

    // create version file with euqal app version and a greater seqan version than the current
    create_file(APP_VERSION_FILENAME, std::string{"2.3.4\n20.5.9"});

    // create timestamp file that dates one day before current to trigger a message (one day = 86400 seconds)
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW - 100401));

    auto [out, err, app_call_succeeded] = simulate_argument_parser(3, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, detail::version_checker::message_seqan3_update);

    EXPECT_TRUE(std::regex_match(read_file(APP_TIMESTAMP_FILENAME), timestamp_regex));

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

// case: the current argument parser has a greater app version than is present in the version file
TEST_F(version_check, greater_app_version)
{
    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_ON};

    // create version file with equal seqan version and a smaller app version than the current
    ASSERT_TRUE(create_file(APP_VERSION_FILENAME, std::string{"1.5.9\n" + SEQAN_VERSION}));

    // create timestamp file that dates one day before current to trigger a message
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW - 100401)); // one day = 86400 seconds

    auto [out, err, app_call_succeeded] = simulate_argument_parser(3, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, detail::version_checker::message_registered_app_update);

    EXPECT_TRUE(std::regex_match(read_file(APP_TIMESTAMP_FILENAME), timestamp_regex));

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, unregistered_app)
{
    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_ON};

    // create version file with equal seqan version and a smaller app version than the current
    ASSERT_TRUE(create_file(APP_VERSION_FILENAME, std::string{ "UNREGISTERED_APP\n" + SEQAN_VERSION}));

    // create timestamp file that dates one day before current to trigger a message
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW - 100401)); // one day = 86400 seconds

    auto [out, err, app_call_succeeded] = simulate_argument_parser(3, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, detail::version_checker::message_unregistered_app);

    EXPECT_TRUE(std::regex_match(read_file(APP_TIMESTAMP_FILENAME), timestamp_regex));

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}
#endif // !defined(NDEBUG)

// case: the current argument parser has a smaller app version than is present in the version file
#if defined(NDEBUG)
TEST_F(version_check, smaller_app_version)
{
    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_ON};

    // create version file with equal seqan version and a greater app version than the current
    ASSERT_TRUE(create_file(APP_VERSION_FILENAME, std::string{"20.5.9\n" + SEQAN_VERSION}));

    // create timestamp file that dates one day before current to trigger a message (one day = 86400 seconds)
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW - 100401));

    auto [out, err, app_call_succeeded] = simulate_argument_parser(3, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, (detail::version_checker{APP_NAME, "2.3.4"}.message_app_update));

    EXPECT_TRUE(std::regex_match(read_file(APP_TIMESTAMP_FILENAME), timestamp_regex));

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, smaller_app_version_custom_url)
{
    char * env{std::getenv("SEQAN3_NO_VERSION_CHECK")};
    if (env != nullptr)
        unsetenv("SEQAN3_NO_VERSION_CHECK");

    const char * argv[3] = {APP_NAME.c_str(), OPTION_VERSION_CHECK, OPTION_ON};

    // create version file with equal seqan version and a greater app version than the current
    ASSERT_TRUE(create_file(APP_VERSION_FILENAME, std::string{"20.5.9\n" + SEQAN_VERSION}));

    // create timestamp file that dates one day before current to trigger a message (one day = 86400 seconds)
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, TIME_NOW - 100401));

    argument_parser parser{APP_NAME, 3, argv};
    parser.info.version = "2.3.4";
    parser.info.url = "https//foo.de";

    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();
    EXPECT_NO_THROW(parser.parse());
    std::string out = testing::internal::GetCapturedStdout();
    std::string err = testing::internal::GetCapturedStderr();

    // call future.get() to artificially wait for the thread to finish and avoid
    // any interference with following tests
    if (parser.version_check_future.valid())
        parser.version_check_future.get();

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, (detail::version_checker{APP_NAME, parser.info.version, parser.info.url}.message_app_update));

    EXPECT_TRUE(std::regex_match(read_file(APP_TIMESTAMP_FILENAME), timestamp_regex));

    if (env != nullptr)
        setenv("SEQAN3_NO_VERSION_CHECK", env, 1);

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}
#endif // defined(NDEBUG)

TEST_F(version_check, user_specified_never)
{
    const char * argv[2] = {APP_NAME.c_str(), "-f"}; // no explicit version check option

    // create timestamp files
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, "NEVER"));

    auto [out, err, app_call_succeeded] = simulate_argument_parser(2, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    EXPECT_FALSE(std::filesystem::exists(APP_VERSION_FILENAME));
    EXPECT_EQ(read_file(APP_TIMESTAMP_FILENAME), "NEVER"); // should not be modified

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, user_specified_always)
{
    const char * argv[2] = {APP_NAME.c_str(), "-f"}; // no explicit version check option

    // create timestamp files
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, "ALWAYS"));

    auto [out, err, app_call_succeeded] = simulate_argument_parser(2, argv);

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    if (app_call_succeeded)
    {
        EXPECT_TRUE(std::filesystem::exists(APP_VERSION_FILENAME));
    }
    else
    {
        std::cout << "App call did not succeed (server offline?) and could thus not be tested." << std::endl;
    }

    EXPECT_EQ(read_file(APP_TIMESTAMP_FILENAME), "ALWAYS"); // should not be modified
    EXPECT_TRUE(remove_files_from_path()); // clear files again

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}

TEST_F(version_check, wrong_version_string)
{
    const char * argv[2] = {APP_NAME.c_str(), "-f"}; // no explicit version check option

    // create a corrupted version file. Nothing should be printed, it is just ignored
    ASSERT_TRUE(create_file(APP_VERSION_FILENAME, std::string{"20.wrong.9\nalso.wrong.4"}));
    ASSERT_TRUE(create_file(APP_TIMESTAMP_FILENAME, "ALWAYS"));

    auto [out, err, app_call_succeeded] = simulate_argument_parser(2, argv);
    (void) app_call_succeeded;

    EXPECT_EQ(out, "");
    EXPECT_EQ(err, "");

    EXPECT_TRUE(remove_files_from_path()); // clear files again
}
