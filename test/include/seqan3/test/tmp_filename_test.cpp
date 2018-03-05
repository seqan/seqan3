// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <gtest/gtest.h>
#include <fstream>
#include "tmp_filename.hpp"

using namespace seqan3::test;

// aggregate initialization
TEST(tmp_filename_aggr, aggr)
{
    tmp_file_name t1{"aggr_test"};
    tmp_file_name t2("aggr_test");
    EXPECT_NE(t1.get_path(), t2.get_path());
    EXPECT_TRUE(fs::exists(t1.get_path().parent_path()));
    EXPECT_TRUE(fs::exists(t2.get_path().parent_path()));
    EXPECT_EQ(fs::temp_directory_path(), t1.get_path().parent_path().parent_path()/="/");
    EXPECT_EQ(fs::temp_directory_path(), t2.get_path().parent_path().parent_path()/="/");
}

// nullptr as filename
TEST(tmp_filename_nullptr, null_ptr)
{
    EXPECT_THROW(tmp_file_name t1{nullptr}, fs::filesystem_error);
    EXPECT_THROW(tmp_file_name t1(nullptr), fs::filesystem_error);
}

// move construction
TEST(tmp_filename_mv_ctr, mv_ctr)
{
    tmp_file_name t1{"mv_ctr_test"};
    tmp_file_name t2{"mv_ctr_test"};
    tmp_file_name t3{std::move(t2)};
    EXPECT_NE(t1.get_path(), t3.get_path());
    tmp_file_name t4(std::move(t1));
    EXPECT_NE(t3.get_path(), t4.get_path());
}

// move assignment
TEST(tmp_filename_mv_assign, mv_assign)
{
    tmp_file_name t1{"mv_ctr_test"};
    tmp_file_name t2{"mv_ctr_test"};
    tmp_file_name t3 = std::move(t2);
    EXPECT_NE(t1.get_path(), t3.get_path());
}

// destructor
TEST(tmp_filename_dtr, dtr)
{
    auto t1 = std::make_unique<tmp_file_name>("delete_test");
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
