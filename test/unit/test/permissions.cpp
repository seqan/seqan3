// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/file_access.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/tmp_filename.hpp>

TEST(read_access, granted)
{
    seqan3::test::tmp_filename file{"seqan3_test_access_read_file_granted"};
    std::filesystem::path path{file.get_path()};
    {
        std::ofstream str{path};
    }
    EXPECT_TRUE(seqan3::test::read_access(path));
}

TEST(read_access, revoked)
{
    seqan3::test::tmp_filename file{"seqan3_test_access_read_file_revoked"};
    std::filesystem::path path{file.get_path()};
    {
        std::ofstream str{path};
    }
    std::filesystem::permissions(path,
                                 std::filesystem::perms::owner_read,
                                 std::filesystem::perm_options::remove);

    auto check_readability = [&] ()
    {
        std::fstream stream;
        stream.open(path, std::ios::in);
        return !stream.fail();
    };

    EXPECT_EQ(seqan3::test::read_access(path), check_readability());
}

TEST(write_access, granted_file)
{
    seqan3::test::tmp_filename file{"seqan3_test_access_write_file_granted"};
    EXPECT_TRUE(seqan3::test::write_access(file.get_path()));
}

TEST(write_access, granted_directory)
{
    seqan3::test::tmp_directory directory{};
    EXPECT_TRUE(seqan3::test::write_access(directory.path()));
}

TEST(write_access, revoked_file)
{
    seqan3::test::tmp_filename file{"seqan3_test_access_write_file_revoked"};
    std::filesystem::path path{file.get_path()};
    {
        std::ofstream str{path};
    }
    std::filesystem::permissions(path,
                                 std::filesystem::perms::owner_write,
                                 std::filesystem::perm_options::remove);

    auto check_writeability = [&] ()
    {
        std::fstream stream;
        stream.open(path, std::ios::out);
        return !stream.fail();
    };

    EXPECT_EQ(seqan3::test::write_access(path), check_writeability());
}

TEST(write_access, revoked_directory)
{
    seqan3::test::tmp_directory directory{};
    std::filesystem::path directory_path{directory.path()};
    std::filesystem::permissions(directory.path(),
                                 std::filesystem::perms::owner_write,
                                 std::filesystem::perm_options::remove);

    std::filesystem::path file_path{directory.path()};
    file_path /= "seqan3_test_write_access_check_writeability";

    auto check_writeability = [&] ()
    {
        std::fstream stream;
        stream.open(file_path, std::ios::out);
        return !stream.fail();
    };

    EXPECT_EQ(seqan3::test::write_access(directory.path()), check_writeability());

    if (std::filesystem::exists(file_path))
        std::filesystem::remove(file_path);
}
