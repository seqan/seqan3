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
#include <sstream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;
using namespace seqan3::literal;

// default/zero construction
TEST(dna4, ctr)
{
    dna4 t1;
}

// zero initialization
TEST(dna4, zro)
{
    dna4 t0{};
    EXPECT_EQ(t0, dna4::A);
}

// copy construction
TEST(dna4, cp_ctr)
{
    dna4 t1{dna4::C};
    dna4 t2{t1};
    dna4 t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(dna4, mv_ctr)
{
    dna4 t0{dna4::C};
    dna4 t1{dna4::C};
    dna4 t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    dna4 t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(dna4, cp_assgn)
{
    dna4 t1{dna4::C};
    dna4 t2;
    dna4 t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(dna4, mv_assgn)
{
    dna4 t0{dna4::C};
    dna4 t1{dna4::C};
    dna4 t2;
    dna4 t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(dna4, swap)
{
    dna4 t0{dna4::C};
    dna4 t1{dna4::C};
    dna4 t2{};
    dna4 t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// comparison
TEST(dna4, cmp)
{
    dna4 t0{dna4::A};
    dna4 t1{dna4::C};
    dna4 t2{dna4::G};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

TEST(dna4, to_char_member)
{
    EXPECT_EQ(dna4::A.to_char(), 'A');
    EXPECT_EQ(dna4::C.to_char(), 'C');
    EXPECT_EQ(dna4::G.to_char(), 'G');
    EXPECT_EQ(dna4::T.to_char(), 'T');
    EXPECT_EQ(dna4::U.to_char(), 'T');
    EXPECT_EQ(dna4::UNKNOWN.to_char(), 'A');
}

TEST(dna4, to_char_free)
{
    EXPECT_EQ(to_char(dna4::A), 'A');
    EXPECT_EQ(to_char(dna4::C), 'C');
    EXPECT_EQ(to_char(dna4::G), 'G');
    EXPECT_EQ(to_char(dna4::T), 'T');
    EXPECT_EQ(to_char(dna4::U), 'T');
    EXPECT_EQ(to_char(dna4::UNKNOWN), 'A');
}

TEST(dna4, to_rank_member)
{
    EXPECT_EQ(dna4::A.to_rank(), 0);
    EXPECT_EQ(dna4::C.to_rank(), 1);
    EXPECT_EQ(dna4::G.to_rank(), 2);
    EXPECT_EQ(dna4::T.to_rank(), 3);
    EXPECT_EQ(dna4::U.to_rank(), 3);
    EXPECT_EQ(dna4::UNKNOWN.to_rank(), 0);
}

TEST(dna4, to_rank_free)
{
    EXPECT_EQ(to_rank(dna4::A), 0);
    EXPECT_EQ(to_rank(dna4::C), 1);
    EXPECT_EQ(to_rank(dna4::G), 2);
    EXPECT_EQ(to_rank(dna4::T), 3);
    EXPECT_EQ(to_rank(dna4::U), 3);
    EXPECT_EQ(to_rank(dna4::UNKNOWN), 0);
}

TEST(dna4, stream_operator)
{
    std::stringstream ss;
    ss << dna4::A << dna4::C << dna4::G << dna4::T << dna4::U << dna4::UNKNOWN;
    EXPECT_EQ(ss.str(), "ACGTTA");
}

TEST(dna4, assign_char_member)
{
    dna4 t0;
    t0.assign_char('A');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    t0.assign_char('C');
    EXPECT_EQ(t0, dna4::C);
    t0.assign_char('G');
    EXPECT_EQ(t0, dna4::G);
    t0.assign_char('T');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);
    t0.assign_char('U');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);

    t0.assign_char('a');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    t0.assign_char('c');
    EXPECT_EQ(t0, dna4::C);
    t0.assign_char('g');
    EXPECT_EQ(t0, dna4::G);
    t0.assign_char('t');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);
    t0.assign_char('u');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);

    t0.assign_char('z');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    t0.assign_char('H');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    t0.assign_char('*');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);

    static_assert(std::is_same_v<decltype(t0.assign_char('C')), dna4 &>);
    EXPECT_EQ(t0.assign_char('C'), dna4::C);
}

TEST(dna4, assign_char_free)
{
    dna4 t0;
    assign_char(t0, 'A');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    assign_char(t0, 'C');
    EXPECT_EQ(t0, dna4::C);
    assign_char(t0, 'G');
    EXPECT_EQ(t0, dna4::G);
    assign_char(t0, 'T');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);
    assign_char(t0, 'U');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);

    assign_char(t0, 'a');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    assign_char(t0, 'c');
    EXPECT_EQ(t0, dna4::C);
    assign_char(t0, 'g');
    EXPECT_EQ(t0, dna4::G);
    assign_char(t0, 't');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);
    assign_char(t0, 'u');
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);

    assign_char(t0, 'z');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    assign_char(t0, 'H');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    assign_char(t0, '*');
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);

    static_assert(std::is_same_v<decltype(assign_char(t0, 'C')), dna4 &>);
    EXPECT_EQ(assign_char(t0, 'C'), dna4::C);
}

TEST(dna4, assign_rank_member)
{
    dna4 t0;
    t0.assign_rank(0);
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    t0.assign_rank(1);
    EXPECT_EQ(t0, dna4::C);
    t0.assign_rank(2);
    EXPECT_EQ(t0, dna4::G);
    t0.assign_rank(3);
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);

    static_assert(std::is_same_v<decltype(t0.assign_rank(2)), dna4 &>);
    EXPECT_EQ(t0.assign_rank(1), dna4::C);
}

TEST(dna4, assign_rank_free)
{
    dna4 t0;
    assign_rank(t0, 0);
    EXPECT_EQ(t0, dna4::A);
    EXPECT_EQ(t0, dna4::UNKNOWN);
    assign_rank(t0, 1);
    EXPECT_EQ(t0, dna4::C);
    assign_rank(t0, 2);
    EXPECT_EQ(t0, dna4::G);
    assign_rank(t0, 3);
    EXPECT_EQ(t0, dna4::T);
    EXPECT_EQ(t0, dna4::U);

    static_assert(std::is_same_v<decltype(assign_rank(t0, 2)), dna4 &>);
    EXPECT_EQ(assign_rank(t0, 1), dna4::C);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(dna4_literals, vector)
{
    dna4_vector v;
    v.resize(5, dna4::A);
    EXPECT_EQ(v, "AAAAA"_dna4);

    std::vector<dna4> w{dna4::A, dna4::C, dna4::G, dna4::T, dna4::U, dna4::UNKNOWN};
    EXPECT_EQ(w, "ACGTTA"_dna4);
}

TEST(dna4_literals, basic_string)
{
    dna4_string v;
    v.resize(5, dna4::A);
    EXPECT_EQ(v, "AAAAA"_dna4s);

    std::basic_string<dna4, std::char_traits<dna4>> w{dna4::A, dna4::C, dna4::G, dna4::T, dna4::U, dna4::UNKNOWN};
    EXPECT_EQ(w, "ACGTTA"_dna4s);
}
