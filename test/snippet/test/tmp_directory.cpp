// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/test/tmp_directory.hpp>

TEST(snippet_tmp_directory, tmp_directory_)
{
    // create a directory folder
    seqan3::test::tmp_directory tmp{};

    // Some function that should creates temporary files and removes them again
    {
        std::ofstream ofs{tmp.path() / "somefile.txt"};
        ofs << "Hello World!";
        ofs.close();

        std::filesystem::remove(tmp.path() / "somefile.txt");
    }

    // check that everything was cleanup properly
    EXPECT_TRUE(tmp.empty());
}
