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

#include <iostream>
#include <string>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

using namespace seqan3;

TEST(debug_stream, basic)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << 'a';
    o.flush();
    EXPECT_EQ(o.str(), "a");

    my_stream << "AGA";
    o.flush();
    EXPECT_EQ(o.str(), "aAGA");

    my_stream << 42;
    o.flush();
    EXPECT_EQ(o.str(), "aAGA42");

    int const i = 7;
    my_stream << i;
    o.flush();
    EXPECT_EQ(o.str(), "aAGA427");
}

TEST(debug_stream, range)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    std::vector<int> vec{1, 4, 5, 7, 32, 321};
    my_stream << vec;
    o.flush();
    EXPECT_EQ(o.str(), "[1,4,5,7,32,321]");

    std::vector<std::vector<int>> const vec2 = {{1, 2, 33}, {22,11}};
    my_stream << vec2;
    o.flush();
    EXPECT_EQ(o.str(), "[1,4,5,7,32,321][[1,2,33],[22,11]]");
}

TEST(debug_stream, alphabet)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << 'A'_dna4;
    o.flush();
    EXPECT_EQ(o.str(), "A");

    dna5 d = 'N'_dna5;
    my_stream << d;
    o.flush();
    EXPECT_EQ(o.str(), "AN");

    dna5 const d2 = 'N'_dna5;
    my_stream << d2;
    o.flush();
    EXPECT_EQ(o.str(), "ANN");
}

TEST(debug_stream, range_of_alphabet)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    // temporary
    my_stream << "AGGATAC"_dna5;
    o.flush();
    EXPECT_EQ(o.str(), "AGGATAC");

    // lvalue ref
    auto d = "AGGATAC"_dna5;
    my_stream << d;
    o.flush();
    EXPECT_EQ(o.str(), "AGGATACAGGATAC");

    // lvalue const ref
    auto const d2 = "AGGATAC"_dna5;
    my_stream << d2;
    o.flush();
    EXPECT_EQ(o.str(), "AGGATACAGGATACAGGATAC");

    concatenated_sequences<bitcompressed_vector<dna5>> const vec2 = {"ACGT"_dna5, "GAGGA"_dna5};
    my_stream << vec2;
    o.flush();
    EXPECT_EQ(o.str(), "AGGATACAGGATACAGGATAC[ACGT,GAGGA]");
}

TEST(debug_stream, std_endl)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    my_stream << "foo" << std::endl << "bar";
    o.flush();
    EXPECT_EQ(o.str(), "foo\nbar");
}

TEST(debug_stream, path)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    filesystem::path p{"my/path/my_file.txt"};

    my_stream << p;
    o.flush();
    EXPECT_EQ(o.str(), "\"my/path/my_file.txt\"");
}

TEST(debug_stream, tuple)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    std::tuple<size_t, std::string> t0{32, "dummy"};
    my_stream << t0;
    o.flush();
    EXPECT_EQ(o.str(), "(32,dummy)");

    std::tuple<size_t> t1{32};
    my_stream << t1;
    o.flush();
    EXPECT_EQ(o.str(), "(32,dummy)(32)");

    std::tuple<size_t, std::pair<size_t, size_t>> t2{2, {3,2}};
    my_stream << t2;
    o.flush();
    EXPECT_EQ(o.str(), "(32,dummy)(32)(2,(3,2))");
}
