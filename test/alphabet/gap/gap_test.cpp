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
// Author: David Heller <david.heller@fu-berlin.de>
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gap.hpp>

using namespace seqan3;

TEST(gap_test, test_alphabet_concept)
{
    EXPECT_TRUE(alphabet_concept<gap>);
}

TEST(gap_test, test_default_initialization)
{
    EXPECT_EQ(gap::GAP.to_char(), '-');
}

TEST(gap_test, test_static_cast)
{
    EXPECT_EQ(static_cast<char>(gap{}), '-');
}

TEST(gap_test, test_relations)
{
    EXPECT_EQ(gap{}, gap{});
    EXPECT_LE(gap{}, gap{});
    EXPECT_GE(gap{}, gap{});
}

TEST(gap_test, test_stream_operator)
{
    std::stringstream ss;
    ss << gap{} << gap{} << gap{};
    EXPECT_EQ(ss.str(), "---");
}

TEST( gap_test, test_assign_char)
{
    EXPECT_EQ(gap{}.assign_char('-'), gap{});
    EXPECT_EQ(gap{}.assign_char('x'), gap{});
}

TEST(gap_test, test_to_rank)
{
    EXPECT_EQ(gap{}.to_rank(), 0);
}

TEST(gap_test, test_assign_rank)
{
    EXPECT_EQ(gap{}.assign_rank(0), gap{});
    EXPECT_EQ(gap{}.assign_rank(13), gap{});
}
