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

#include <seqan3/alphabet/nucleotide/rna4.hpp>

using namespace seqan3;
using namespace seqan3::literal;

// default/zero construction
TEST(rna4, ctr)
{
    rna4 t1;
}

// zero initialization
TEST(rna4, zro)
{
    rna4 t0{};
    EXPECT_EQ(t0, rna4::A);
}

// copy construction
TEST(rna4, cp_ctr)
{
    rna4 t1{rna4::C};
    rna4 t2{t1};
    rna4 t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(rna4, mv_ctr)
{
    rna4 t0{rna4::C};
    rna4 t1{rna4::C};
    rna4 t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    rna4 t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(rna4, cp_assgn)
{
    rna4 t1{rna4::C};
    rna4 t2;
    rna4 t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(rna4, mv_assgn)
{
    rna4 t0{rna4::C};
    rna4 t1{rna4::C};
    rna4 t2;
    rna4 t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(rna4, swap)
{
    rna4 t0{rna4::C};
    rna4 t1{rna4::C};
    rna4 t2{};
    rna4 t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// comparison
TEST(rna4, cmp)
{
    rna4 t0{rna4::A};
    rna4 t1{rna4::C};
    rna4 t2{rna4::G};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

TEST(rna4, to_char_member)
{
    EXPECT_EQ(rna4::A.to_char(), 'A');
    EXPECT_EQ(rna4::C.to_char(), 'C');
    EXPECT_EQ(rna4::G.to_char(), 'G');
    EXPECT_EQ(rna4::T.to_char(), 'U');
    EXPECT_EQ(rna4::U.to_char(), 'U');
    EXPECT_EQ(rna4::UNKNOWN.to_char(), 'A');
}

TEST(rna4, to_char_free)
{
    EXPECT_EQ(to_char(rna4::A), 'A');
    EXPECT_EQ(to_char(rna4::C), 'C');
    EXPECT_EQ(to_char(rna4::G), 'G');
    EXPECT_EQ(to_char(rna4::T), 'U');
    EXPECT_EQ(to_char(rna4::U), 'U');
    EXPECT_EQ(to_char(rna4::UNKNOWN), 'A');
}

TEST(rna4, to_rank_member)
{
    EXPECT_EQ(rna4::A.to_rank(), 0);
    EXPECT_EQ(rna4::C.to_rank(), 1);
    EXPECT_EQ(rna4::G.to_rank(), 2);
    EXPECT_EQ(rna4::T.to_rank(), 3);
    EXPECT_EQ(rna4::U.to_rank(), 3);
    EXPECT_EQ(rna4::UNKNOWN.to_rank(), 0);
}

TEST(rna4, to_rank_free)
{
    EXPECT_EQ(to_rank(rna4::A), 0);
    EXPECT_EQ(to_rank(rna4::C), 1);
    EXPECT_EQ(to_rank(rna4::G), 2);
    EXPECT_EQ(to_rank(rna4::T), 3);
    EXPECT_EQ(to_rank(rna4::U), 3);
    EXPECT_EQ(to_rank(rna4::UNKNOWN), 0);
}

TEST(rna4, stream_operator)
{
    std::stringstream ss;
    ss << rna4::A << rna4::C << rna4::G << rna4::T << rna4::U << rna4::UNKNOWN;
    EXPECT_EQ(ss.str(), "ACGUUA");
}

TEST(rna4, assign_char_member)
{
    rna4 t0;
    t0.assign_char('A');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    t0.assign_char('C');
    EXPECT_EQ(t0, rna4::C);
    t0.assign_char('G');
    EXPECT_EQ(t0, rna4::G);
    t0.assign_char('T');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);
    t0.assign_char('U');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);

    t0.assign_char('a');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    t0.assign_char('c');
    EXPECT_EQ(t0, rna4::C);
    t0.assign_char('g');
    EXPECT_EQ(t0, rna4::G);
    t0.assign_char('t');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);
    t0.assign_char('u');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);

    t0.assign_char('z');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    t0.assign_char('H');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    t0.assign_char('*');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);

    static_assert(std::is_same_v<decltype(t0.assign_char('C')), rna4 &>);
    EXPECT_EQ(t0.assign_char('C'), rna4::C);
}

TEST(rna4, assign_char_free)
{
    rna4 t0;
    assign_char(t0, 'A');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    assign_char(t0, 'C');
    EXPECT_EQ(t0, rna4::C);
    assign_char(t0, 'G');
    EXPECT_EQ(t0, rna4::G);
    assign_char(t0, 'T');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);
    assign_char(t0, 'U');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);

    assign_char(t0, 'a');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    assign_char(t0, 'c');
    EXPECT_EQ(t0, rna4::C);
    assign_char(t0, 'g');
    EXPECT_EQ(t0, rna4::G);
    assign_char(t0, 't');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);
    assign_char(t0, 'u');
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);

    assign_char(t0, 'z');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    assign_char(t0, 'H');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    assign_char(t0, '*');
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);

    static_assert(std::is_same_v<decltype(assign_char(t0, 'C')), rna4 &>);
    EXPECT_EQ(assign_char(t0, 'C'), rna4::C);
}

TEST(rna4, assign_rank_member)
{
    rna4 t0;
    t0.assign_rank(0);
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    t0.assign_rank(1);
    EXPECT_EQ(t0, rna4::C);
    t0.assign_rank(2);
    EXPECT_EQ(t0, rna4::G);
    t0.assign_rank(3);
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);

    static_assert(std::is_same_v<decltype(t0.assign_rank(2)), rna4 &>);
    EXPECT_EQ(t0.assign_rank(1), rna4::C);
}

TEST(rna4, assign_rank_free)
{
    rna4 t0;
    assign_rank(t0, 0);
    EXPECT_EQ(t0, rna4::A);
    EXPECT_EQ(t0, rna4::UNKNOWN);
    assign_rank(t0, 1);
    EXPECT_EQ(t0, rna4::C);
    assign_rank(t0, 2);
    EXPECT_EQ(t0, rna4::G);
    assign_rank(t0, 3);
    EXPECT_EQ(t0, rna4::T);
    EXPECT_EQ(t0, rna4::U);

    static_assert(std::is_same_v<decltype(assign_rank(t0, 2)), rna4 &>);
    EXPECT_EQ(assign_rank(t0, 1), rna4::C);
}

// ------------------------------------------------------------------
// compatibility to dna4
// ------------------------------------------------------------------

TEST(rna4, to_dna4)
{
    rna4 t0{rna4::C};
    dna4 d0{t0};
    EXPECT_EQ(d0, dna4::C);

    d0 = rna4::U;
    EXPECT_EQ(d0, dna4::T);
}

TEST(rna4, from_dna4)
{
    dna4 d0{dna4::C};
    rna4 t0{d0};
    EXPECT_EQ(t0, rna4::C);

    t0 = dna4::T;
    EXPECT_EQ(t0, rna4::U);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(rna4_literals, vector)
{
    rna4_vector v;
    v.resize(5, rna4::A);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<rna4> w{rna4::A, rna4::C, rna4::G, rna4::T, rna4::U, rna4::UNKNOWN};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

TEST(rna4_literals, basic_string)
{
    rna4_string v;
    v.resize(5, rna4::A);
    EXPECT_EQ(v, "AAAAA"_rna4s);

    std::basic_string<rna4, std::char_traits<rna4>> w{rna4::A, rna4::C, rna4::G, rna4::T, rna4::U, rna4::UNKNOWN};
    EXPECT_EQ(w, "ACGUUA"_rna4s);
}
