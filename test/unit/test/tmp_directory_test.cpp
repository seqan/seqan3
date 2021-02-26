// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2021-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2021-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <fstream>

#include <seqan3/test/tmp_directory.hpp>

namespace fs = std::filesystem;

// check unique creation of paths
TEST(tmp_directory, unique)
{
    seqan3::test::tmp_directory t1;
    seqan3::test::tmp_directory t2;

    // checking they are unique
    EXPECT_NE(t1.get_path(), t2.get_path());

    // checking they got created
    EXPECT_TRUE(fs::exists(t1.get_path()));
    EXPECT_TRUE(fs::exists(t2.get_path()));

    // check if created paths are empty
    EXPECT_TRUE(t1.empty());
    EXPECT_TRUE(t2.empty());


    // checking they are inside of /tmp
    EXPECT_TRUE(fs::equivalent(fs::temp_directory_path(), fs::path{t1.get_path()}.parent_path()));
    EXPECT_TRUE(fs::equivalent(fs::temp_directory_path(), fs::path{t2.get_path()}.parent_path()));
}

// move construction
TEST(tmp_directory_mv_ctr, mv_ctr)
{
    seqan3::test::tmp_directory t1{};
    seqan3::test::tmp_directory t2{};
    seqan3::test::tmp_directory t3{std::move(t2)};

    EXPECT_TRUE(fs::exists(t1.get_path()));
    EXPECT_TRUE(fs::exists(t3.get_path()));
    EXPECT_TRUE(t1.empty());
    EXPECT_TRUE(t3.empty());

    EXPECT_NE(t1.get_path(), t3.get_path());
    seqan3::test::tmp_directory t4(std::move(t1));

    EXPECT_TRUE(fs::exists(t4.get_path()));
    EXPECT_NE(t3.get_path(), t4.get_path());
}

// move assignment
TEST(tmp_directory_mv_assign, mv_assign)
{
    seqan3::test::tmp_directory t1{};
    seqan3::test::tmp_directory t2{};
    seqan3::test::tmp_directory t3;
    t3 = std::move(t2);

    EXPECT_NE(t1.get_path(), t3.get_path());

    EXPECT_TRUE(fs::exists(t1.get_path()));
    EXPECT_TRUE(fs::exists(t3.get_path()));
}

// check dtor does all its cleanups
TEST(tmp_directory_dtor, dtor)
{
    auto t1 = std::make_unique<seqan3::test::tmp_directory>();
    auto path = t1->get_path();

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
        fs::create_directory(path / "somefolder");
        std::ofstream os{path / "somefolder/file2"};
        os << "other data";
    }
    // check that paths are not empty
    EXPECT_FALSE(t1->empty());

    EXPECT_TRUE(fs::exists(path));
    EXPECT_TRUE(fs::exists(path / "file1"));
    EXPECT_TRUE(fs::exists(path / "somefolder/file2"));

    // Should warn about unclean temporary directory
    testing::internal::CaptureStderr();
    t1.reset();
    std::string output = testing::internal::GetCapturedStderr();
    EXPECT_FALSE(output.empty());

    EXPECT_FALSE(fs::exists(path));
    EXPECT_FALSE(fs::exists(path / "file1"));
    EXPECT_FALSE(fs::exists(path / "somefolder/file2"));
}
