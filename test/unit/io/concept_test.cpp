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

#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <seqan3/io/stream/concept.hpp>

using namespace seqan3;

TEST(io, Stream)
{
    EXPECT_FALSE((Stream<std::istream, char>));
    EXPECT_FALSE((Stream<std::ostream, char>));
    EXPECT_TRUE((Stream<std::iostream, char>));
    EXPECT_FALSE((Stream<std::ifstream, char>));
    EXPECT_FALSE((Stream<std::ofstream, char>));
    EXPECT_TRUE((Stream<std::fstream, char>));
    EXPECT_FALSE((Stream<std::istringstream, char>));
    EXPECT_FALSE((Stream<std::ostringstream, char>));
    EXPECT_TRUE((Stream<std::stringstream, char>));
    EXPECT_FALSE((Stream<std::vector<char>, char>));
    EXPECT_FALSE((Stream<std::vector<char>, char>));
    EXPECT_FALSE((Stream<std::vector<char>, char>));
    EXPECT_FALSE((Stream<std::string, char>));
    EXPECT_FALSE((Stream<std::string, char>));
    EXPECT_FALSE((Stream<std::string, char>));
}

TEST(io, Istream)
{
    EXPECT_TRUE((Istream<std::istream, char>));
    EXPECT_FALSE((Istream<std::ostream, char>));
    EXPECT_TRUE((Istream<std::iostream, char>));
    EXPECT_TRUE((Istream<std::ifstream, char>));
    EXPECT_FALSE((Istream<std::ofstream, char>));
    EXPECT_TRUE((Istream<std::fstream, char>));
    EXPECT_TRUE((Istream<std::istringstream, char>));
    EXPECT_FALSE((Istream<std::ostringstream, char>));
    EXPECT_TRUE((Istream<std::stringstream, char>));
    EXPECT_FALSE((Istream<std::vector<char>, char>));
    EXPECT_FALSE((Istream<std::vector<char>, char>));
    EXPECT_FALSE((Istream<std::vector<char>, char>));
    EXPECT_FALSE((Istream<std::string, char>));
    EXPECT_FALSE((Istream<std::string, char>));
    EXPECT_FALSE((Istream<std::string, char>));
}

TEST(io, Ostream)
{
    EXPECT_FALSE((Ostream<std::istream, char>));
    EXPECT_TRUE((Ostream<std::ostream, char>));
    EXPECT_TRUE((Ostream<std::iostream, char>));
    EXPECT_FALSE((Ostream<std::ifstream, char>));
    EXPECT_TRUE((Ostream<std::ofstream, char>));
    EXPECT_TRUE((Ostream<std::fstream, char>));
    EXPECT_FALSE((Ostream<std::istringstream, char>));
    EXPECT_TRUE((Ostream<std::ostringstream, char>));
    EXPECT_TRUE((Ostream<std::stringstream, char>));
    EXPECT_FALSE((Ostream<std::vector<char>, char>));
    EXPECT_FALSE((Ostream<std::vector<char>, char>));
    EXPECT_FALSE((Ostream<std::vector<char>, char>));
    EXPECT_FALSE((Ostream<std::string, char>));
    EXPECT_FALSE((Ostream<std::string, char>));
    EXPECT_FALSE((Ostream<std::string, char>));
}
