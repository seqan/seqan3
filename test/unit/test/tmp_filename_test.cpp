// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

namespace fs = std::filesystem;

// aggregate initialization
TEST(tmp_filename_aggr, aggr)
{
    seqan3::test::tmp_filename t1{"aggr_test"};
    seqan3::test::tmp_filename t2("aggr_test");
    EXPECT_NE(t1.get_path(), t2.get_path());
    EXPECT_TRUE(fs::exists(t1.get_path().parent_path()));
    EXPECT_TRUE(fs::exists(t2.get_path().parent_path()));

    EXPECT_TRUE(fs::equivalent(fs::temp_directory_path(), t1.get_path().parent_path().parent_path()));
    EXPECT_TRUE(fs::equivalent(fs::temp_directory_path(), t2.get_path().parent_path().parent_path()));
}

// nullptr as filename
TEST(tmp_filename_nullptr, null_ptr)
{
    EXPECT_THROW(seqan3::test::tmp_filename t1{nullptr}, fs::filesystem_error);
    EXPECT_THROW(seqan3::test::tmp_filename t1(nullptr), fs::filesystem_error);
}

// move construction
TEST(tmp_filename_mv_ctr, mv_ctr)
{
    seqan3::test::tmp_filename t1{"mv_ctr_test"};
    seqan3::test::tmp_filename t2{"mv_ctr_test"};
    seqan3::test::tmp_filename t3{std::move(t2)};
    EXPECT_NE(t1.get_path(), t3.get_path());
    seqan3::test::tmp_filename t4(std::move(t1));
    EXPECT_NE(t3.get_path(), t4.get_path());
}

// move assignment
TEST(tmp_filename_mv_assign, mv_assign)
{
    seqan3::test::tmp_filename t1{"mv_ctr_test"};
    seqan3::test::tmp_filename t2{"mv_ctr_test"};
    seqan3::test::tmp_filename t3 = std::move(t2);
    EXPECT_NE(t1.get_path(), t3.get_path());
}

// destructor
TEST(tmp_filename_dtr, dtr)
{
    auto t1 = std::make_unique<seqan3::test::tmp_filename>("delete_test");
    auto path = t1->get_path();
    std::ofstream os{path, std::ios::out};
    os << "delete_test";
    os.close();
    EXPECT_TRUE(fs::exists(path));
    EXPECT_TRUE(fs::exists(path.parent_path()));
    t1.reset();
    EXPECT_FALSE(fs::exists(path));
    EXPECT_FALSE(fs::exists(path.parent_path()));
}

// throw if invalid TMPDIR
TEST(tmp_filename_throw, throw)
{
    if (!system("man putenv > /dev/null 2>&1"))
    {
        char str[] = "TMPDIR=/invalid";
        putenv(str);
        EXPECT_THROW(seqan3::test::tmp_filename t1{"throw"}, std::filesystem::filesystem_error);
    }
}
