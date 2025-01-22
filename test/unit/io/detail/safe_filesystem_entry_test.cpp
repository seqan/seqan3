// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/io/detail/safe_filesystem_entry.hpp>
#include <seqan3/test/tmp_directory.hpp>

TEST(safe_filesystem_entry, file)
{
    seqan3::test::tmp_directory tmp;
    auto my_file = tmp.path() / "dummy.txt";
    {
        std::ofstream file{my_file};
        EXPECT_TRUE(std::filesystem::exists(my_file));
        seqan3::detail::safe_filesystem_entry file_guard{my_file};
    }

    EXPECT_FALSE(std::filesystem::exists(my_file));
}

TEST(safe_filesystem_entry, file_already_removed)
{
    seqan3::test::tmp_directory tmp;
    auto my_file = tmp.path() / "dummy.txt";
    {
        EXPECT_FALSE(std::filesystem::exists(my_file));
        seqan3::detail::safe_filesystem_entry file_guard{my_file};
    }

    EXPECT_FALSE(std::filesystem::exists(my_file));
}

TEST(safe_filesystem_entry, directory)
{
    seqan3::test::tmp_directory tmp;
    auto my_dir = tmp.path() / "dummy.txt";
    {
        std::filesystem::create_directory(my_dir);
        EXPECT_TRUE(std::filesystem::exists(my_dir));
        seqan3::detail::safe_filesystem_entry dir_guard{my_dir};
    }

    EXPECT_FALSE(std::filesystem::exists(my_dir));
}

TEST(safe_filesystem_entry, directory_already_removed)
{
    seqan3::test::tmp_directory tmp;
    auto my_dir = tmp.path() / "dummy.txt";
    {
        EXPECT_FALSE(std::filesystem::exists(my_dir));
        seqan3::detail::safe_filesystem_entry dir_guard{my_dir};
    }

    EXPECT_FALSE(std::filesystem::exists(my_dir));
}

TEST(safe_filesystem_entry, remove)
{
    seqan3::test::tmp_directory tmp;
    auto my_file = tmp.path() / "dummy.txt";
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
    seqan3::test::tmp_directory tmp;
    auto my_dir = tmp.path() / "dummy.txt";
    {
        std::filesystem::create_directory(my_dir);
        EXPECT_TRUE(std::filesystem::exists(my_dir));
        seqan3::detail::safe_filesystem_entry dir_guard{my_dir};
        EXPECT_TRUE(dir_guard.remove_all());
    }

    EXPECT_FALSE(std::filesystem::exists(my_dir));
}
