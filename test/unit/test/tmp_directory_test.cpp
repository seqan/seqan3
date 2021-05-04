// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2021-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2021-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <fstream>
#include <gtest/gtest.h>
#include <regex>

#include <seqan3/test/tmp_directory.hpp>

// check unique creation of paths
TEST(tmp_directory, unique)
{
    seqan3::test::tmp_directory t1;
    seqan3::test::tmp_directory t2;

    // checking they are unique
    EXPECT_NE(t1.path(), t2.path());

    // checking they got created
    EXPECT_TRUE(std::filesystem::exists(t1.path()));
    EXPECT_TRUE(std::filesystem::exists(t2.path()));

    // check if created folders are empty
    EXPECT_TRUE(t1.empty());
    EXPECT_TRUE(t2.empty());

    // checking they are inside of /tmp
    EXPECT_TRUE(std::filesystem::equivalent(std::filesystem::temp_directory_path(), std::filesystem::path{t1.path()}.parent_path()));
    EXPECT_TRUE(std::filesystem::equivalent(std::filesystem::temp_directory_path(), std::filesystem::path{t2.path()}.parent_path()));
}

// move construction
TEST(tmp_directory, move_constructible)
{
    seqan3::test::tmp_directory t1{};
    seqan3::test::tmp_directory t2{};
    seqan3::test::tmp_directory t3{std::move(t2)};

    EXPECT_TRUE(std::filesystem::exists(t1.path()));
    EXPECT_TRUE(std::filesystem::exists(t3.path()));
    EXPECT_TRUE(t1.empty());
    EXPECT_TRUE(t3.empty());

    EXPECT_NE(t1.path(), t3.path());
    seqan3::test::tmp_directory t4(std::move(t1));

    EXPECT_TRUE(std::filesystem::exists(t4.path()));
    EXPECT_NE(t3.path(), t4.path());
}

// move assignment
TEST(tmp_directory, move_assignable)
{
    std::filesystem::path p1;
    std::filesystem::path p2;
    std::filesystem::path p3;

    {
        seqan3::test::tmp_directory t1{};
        seqan3::test::tmp_directory t2{};
        seqan3::test::tmp_directory t3;

        p1 = t1.path();
        p2 = t2.path();
        p3 = t3.path();

        t3 = std::move(t2);

        EXPECT_NE(t1.path(), t3.path());

        EXPECT_TRUE(std::filesystem::exists(t1.path()));
        EXPECT_TRUE(std::filesystem::exists(t3.path()));
    }
    // check all temporary directories are cleaned
    EXPECT_FALSE(std::filesystem::exists(p1));
    EXPECT_FALSE(std::filesystem::exists(p2));
    EXPECT_FALSE(std::filesystem::exists(p3));

}

// check destructor does all its cleanups
TEST(tmp_directory, cleanup_on_destruction)
{
    std::filesystem::path path;
    {
        seqan3::test::tmp_directory t1{};
        path = t1.path();

        // create file structure
        // /tmp
        //  + seqan3_test_XXXXXXXX
        //    - file1
        //    + somefolder
        //      - file2
        //
        // create file1
        {
            std::ofstream os{path / "file1", std::ios::out};
            os << "some data";
        }
        // create somefoleder/file2
        {
            std::filesystem::create_directory(path / "somefolder");
            std::ofstream os{path / "somefolder/file2"};
            os << "other data";
        }
        // check that paths are not empty
        EXPECT_FALSE(t1.empty());

        EXPECT_TRUE(std::filesystem::exists(path));
        EXPECT_TRUE(std::filesystem::exists(path / "file1"));
        EXPECT_TRUE(std::filesystem::exists(path / "somefolder/file2"));

        // Should warn about unclean temporary directory
        testing::internal::CaptureStderr();
    }

    std::string output = testing::internal::GetCapturedStderr();
    EXPECT_FALSE(output.empty());

    std::regex re{"temporary directory \"[^ ]*\" has some files that should be deleted\n"
                  "- \"[^ ]*\"\n"
                  "- \"[^ ]*\"\n"
                  "- \"[^ ]*\"\n"};
    EXPECT_TRUE(std::regex_match(output, re)) << "Actual output: " << output;

    EXPECT_FALSE(std::filesystem::exists(path));
    EXPECT_FALSE(std::filesystem::exists(path / "file1"));
    EXPECT_FALSE(std::filesystem::exists(path / "somefolder/file2"));
}

// check destructor warns if someone else deletes the temp directory
TEST(tmp_directory, warn_about_missing_managed_tmp_directory_on_destruction)
{
    std::filesystem::path path;
    {
        seqan3::test::tmp_directory t1{};
        path = t1.path();

        // create file structure
        // /tmp
        //  + seqan3_test_XXXXXXXX

        std::filesystem::remove_all(t1.path());

        // Should warn about unclean temporary directory
        testing::internal::CaptureStderr();
    }

    std::string output = testing::internal::GetCapturedStderr();

    std::regex re{"temporary directory \"[^ ]*\" was deleted externally. This is discouraged program behaviour\n"};
    EXPECT_TRUE(std::regex_match(output, re)) << "Actual output: " << output;
}

// check a unwritable tmp file fails
TEST(tmp_directory_throw, directory_not_writeable)
{
    // create a temporary folder that will mimic the normal tmp folder
    seqan3::test::tmp_directory temporary_tmp_folder;
    setenv("TMPDIR", temporary_tmp_folder.path().c_str(), 1); // name, value, overwrite

    // make temporary_tmp_folder read only
    std::filesystem::permissions(temporary_tmp_folder.path(),
                                 std::filesystem::perms::owner_write,
                                 std::filesystem::perm_options::remove);

    // The actual test
    EXPECT_THROW(seqan3::test::tmp_directory{}, std::filesystem::filesystem_error);

    // give temporary_tmp_folder write permissions back
    std::filesystem::permissions(temporary_tmp_folder.path(),
                                 std::filesystem::perms::owner_write,
                                 std::filesystem::perm_options::add);

}
