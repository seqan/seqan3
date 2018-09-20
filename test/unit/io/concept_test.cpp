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

<<<<<<< HEAD
#include <seqan3/io/concept.hpp>
=======
#include <seqan3/io/stream/concept.hpp>
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

using namespace seqan3;

TEST(io, stream_concept)
{
    EXPECT_FALSE((stream_concept<std::istream, char>));
    EXPECT_FALSE((stream_concept<std::ostream, char>));
    EXPECT_TRUE((stream_concept<std::iostream, char>));
    EXPECT_FALSE((stream_concept<std::ifstream, char>));
    EXPECT_FALSE((stream_concept<std::ofstream, char>));
    EXPECT_TRUE((stream_concept<std::fstream, char>));
    EXPECT_FALSE((stream_concept<std::istringstream, char>));
    EXPECT_FALSE((stream_concept<std::ostringstream, char>));
    EXPECT_TRUE((stream_concept<std::stringstream, char>));
    EXPECT_FALSE((stream_concept<std::vector<char>, char>));
    EXPECT_FALSE((stream_concept<std::vector<char>, char>));
    EXPECT_FALSE((stream_concept<std::vector<char>, char>));
    EXPECT_FALSE((stream_concept<std::string, char>));
    EXPECT_FALSE((stream_concept<std::string, char>));
    EXPECT_FALSE((stream_concept<std::string, char>));
}

TEST(io, istream_concept)
{
    EXPECT_TRUE((istream_concept<std::istream, char>));
    EXPECT_FALSE((istream_concept<std::ostream, char>));
    EXPECT_TRUE((istream_concept<std::iostream, char>));
    EXPECT_TRUE((istream_concept<std::ifstream, char>));
    EXPECT_FALSE((istream_concept<std::ofstream, char>));
    EXPECT_TRUE((istream_concept<std::fstream, char>));
    EXPECT_TRUE((istream_concept<std::istringstream, char>));
    EXPECT_FALSE((istream_concept<std::ostringstream, char>));
    EXPECT_TRUE((istream_concept<std::stringstream, char>));
    EXPECT_FALSE((istream_concept<std::vector<char>, char>));
    EXPECT_FALSE((istream_concept<std::vector<char>, char>));
    EXPECT_FALSE((istream_concept<std::vector<char>, char>));
    EXPECT_FALSE((istream_concept<std::string, char>));
    EXPECT_FALSE((istream_concept<std::string, char>));
    EXPECT_FALSE((istream_concept<std::string, char>));
}

TEST(io, ostream_concept)
{
    EXPECT_FALSE((ostream_concept<std::istream, char>));
    EXPECT_TRUE((ostream_concept<std::ostream, char>));
    EXPECT_TRUE((ostream_concept<std::iostream, char>));
    EXPECT_FALSE((ostream_concept<std::ifstream, char>));
    EXPECT_TRUE((ostream_concept<std::ofstream, char>));
    EXPECT_TRUE((ostream_concept<std::fstream, char>));
    EXPECT_FALSE((ostream_concept<std::istringstream, char>));
    EXPECT_TRUE((ostream_concept<std::ostringstream, char>));
    EXPECT_TRUE((ostream_concept<std::stringstream, char>));
    EXPECT_FALSE((ostream_concept<std::vector<char>, char>));
    EXPECT_FALSE((ostream_concept<std::vector<char>, char>));
    EXPECT_FALSE((ostream_concept<std::vector<char>, char>));
    EXPECT_FALSE((ostream_concept<std::string, char>));
    EXPECT_FALSE((ostream_concept<std::string, char>));
    EXPECT_FALSE((ostream_concept<std::string, char>));
}
