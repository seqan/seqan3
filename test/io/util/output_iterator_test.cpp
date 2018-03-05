// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

#include <vector>
#include <sstream>

#include <seqan3/io/util/output_iterator.hpp>

using namespace seqan3;

TEST(output_iterator, vector)
{
    std::string input{"12345 6789"};
    std::vector<char> vec;
    auto it = single_pass_output_iterator(vec);
    for (auto val : input)
    {
        *it = val;
        ++it;
    }
    EXPECT_EQ(vec[0], '1');
    EXPECT_EQ(vec[1], '2');
    EXPECT_EQ(vec[2], '3');
    EXPECT_EQ(vec[3], '4');
    EXPECT_EQ(vec[4], '5');
    EXPECT_EQ(vec[5], ' ');
    EXPECT_EQ(vec[6], '6');
    EXPECT_EQ(vec[7], '7');
    EXPECT_EQ(vec[8], '8');
    EXPECT_EQ(vec[9], '9');
}

TEST(output_iterator, ostream)
{
    std::string input{"12345 6789"};
    std::ostringstream stream;
    auto it = single_pass_output_iterator(stream);
    for (auto val : input)
    {
        *it = val;
        ++it;
    }

    EXPECT_EQ(stream.str(), "12345 6789");
}
