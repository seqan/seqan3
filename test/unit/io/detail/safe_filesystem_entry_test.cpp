// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/io/detail/safe_filesystem_entry.hpp>
#include <seqan3/test/tmp_filename.hpp>

TEST(safe_filesystem_entry, file)
{
    seqan3::test::tmp_filename tmp_file{"dummy.txt"};
    std::filesystem::path my_file = tmp_file.get_path();
    {
        std::ofstream file{my_file};
        EXPECT_TRUE(std::filesystem::exists(my_file));
        seqan3::detail::safe_filesystem_entry file_guard{my_file};
    }

    EXPECT_FALSE(std::filesystem::exists(my_file));
}

TEST(safe_filesystem_entry, file_already_removed)
{
    seqan3::test::tmp_filename tmp_file{"dummy.txt"};
    std::filesystem::path my_file = tmp_file.get_path();
    {
        EXPECT_FALSE(std::filesystem::exists(my_file));
        seqan3::detail::safe_filesystem_entry file_guard{my_file};
    }

    EXPECT_FALSE(std::filesystem::exists(my_file));
}

TEST(safe_filesystem_entry, directory)
{
    seqan3::test::tmp_filename tmp_file{"dummy.txt"};
    std::filesystem::path my_dir = tmp_file.get_path();
    {
        std::filesystem::create_directory(my_dir);
        EXPECT_TRUE(std::filesystem::exists(my_dir));
        seqan3::detail::safe_filesystem_entry dir_guard{my_dir};
    }

    EXPECT_FALSE(std::filesystem::exists(my_dir));
}

TEST(safe_filesystem_entry, directory_already_removed)
{
    seqan3::test::tmp_filename tmp_file{"dummy.txt"};
    std::filesystem::path my_dir = tmp_file.get_path();
    {
        EXPECT_FALSE(std::filesystem::exists(my_dir));
        seqan3::detail::safe_filesystem_entry dir_guard{my_dir};
    }

    EXPECT_FALSE(std::filesystem::exists(my_dir));
}

TEST(safe_filesystem_entry, remove)
{
    seqan3::test::tmp_filename tmp_file{"dummy.txt"};
    std::filesystem::path my_file = tmp_file.get_path();
    {
        std::ofstream file{my_file};
        EXPECT_TRUE(std::filesystem::exists(my_file));
        seqan3::detail::safe_filesystem_entry file_guard{my_file};
        EXPECT_TRUE(file_guard.remove());
    }

    EXPECT_FALSE(std::filesystem::exists(my_file));
}

TEST(safe_filesystem_entry, remove_all)
{
    seqan3::test::tmp_filename tmp_file{"dummy.txt"};
    std::filesystem::path my_dir = tmp_file.get_path();
    {
        std::filesystem::create_directory(my_dir);
        EXPECT_TRUE(std::filesystem::exists(my_dir));
        seqan3::detail::safe_filesystem_entry dir_guard{my_dir};
        EXPECT_TRUE(dir_guard.remove_all());
    }

    EXPECT_FALSE(std::filesystem::exists(my_dir));
}
