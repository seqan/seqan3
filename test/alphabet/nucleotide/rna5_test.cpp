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

#include <seqan3/alphabet/nucleotide/rna5.hpp>

using namespace seqan3;
using namespace seqan3::literal;

// default/zero construction
TEST(rna5, ctr)
{
    rna5 t1;
}

// zero initialization
TEST(rna5, zro)
{
    rna5 t0{};
    EXPECT_EQ(t0, rna5::A);
}

// copy construction
TEST(rna5, cp_ctr)
{
    rna5 t1{rna5::C};
    rna5 t2{t1};
    rna5 t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(rna5, mv_ctr)
{
    rna5 t0{rna5::C};
    rna5 t1{rna5::C};
    rna5 t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    rna5 t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(rna5, cp_assgn)
{
    rna5 t1{rna5::C};
    rna5 t2;
    rna5 t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(rna5, mv_assgn)
{
    rna5 t0{rna5::C};
    rna5 t1{rna5::C};
    rna5 t2;
    rna5 t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(rna5, swap)
{
    rna5 t0{rna5::C};
    rna5 t1{rna5::C};
    rna5 t2{};
    rna5 t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// comparison
TEST(rna5, cmp)
{
    rna5 t0{rna5::A};
    rna5 t1{rna5::C};
    rna5 t2{rna5::G};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

TEST(rna5, to_char_member)
{
    EXPECT_EQ(rna5::A.to_char(), 'A');
    EXPECT_EQ(rna5::C.to_char(), 'C');
    EXPECT_EQ(rna5::G.to_char(), 'G');
    EXPECT_EQ(rna5::T.to_char(), 'U');
    EXPECT_EQ(rna5::U.to_char(), 'U');
    EXPECT_EQ(rna5::N.to_char(), 'N');
    EXPECT_EQ(rna5::UNKNOWN.to_char(), 'N');
}

TEST(rna5, to_char_free)
{
    EXPECT_EQ(to_char(rna5::A), 'A');
    EXPECT_EQ(to_char(rna5::C), 'C');
    EXPECT_EQ(to_char(rna5::G), 'G');
    EXPECT_EQ(to_char(rna5::T), 'U');
    EXPECT_EQ(to_char(rna5::U), 'U');
    EXPECT_EQ(to_char(rna5::N), 'N');
    EXPECT_EQ(to_char(rna5::UNKNOWN), 'N');
}

TEST(rna5, to_rank_member)
{
    EXPECT_EQ(rna5::A.to_rank(), 0);
    EXPECT_EQ(rna5::C.to_rank(), 1);
    EXPECT_EQ(rna5::G.to_rank(), 2);
    EXPECT_EQ(rna5::T.to_rank(), 3);
    EXPECT_EQ(rna5::U.to_rank(), 3);
    EXPECT_EQ(rna5::N.to_rank(), 4);
    EXPECT_EQ(rna5::UNKNOWN.to_rank(), 4);
}

TEST(rna5, to_rank_free)
{
    EXPECT_EQ(to_rank(rna5::A), 0);
    EXPECT_EQ(to_rank(rna5::C), 1);
    EXPECT_EQ(to_rank(rna5::G), 2);
    EXPECT_EQ(to_rank(rna5::T), 3);
    EXPECT_EQ(to_rank(rna5::U), 3);
    EXPECT_EQ(to_rank(rna5::N), 4);
    EXPECT_EQ(to_rank(rna5::UNKNOWN), 4);
}

TEST(rna5, stream_operator)
{
    std::stringstream ss;
    ss << rna5::A << rna5::C << rna5::G << rna5::T << rna5::U << rna5::N << rna5::UNKNOWN;
    EXPECT_EQ(ss.str(), "ACGUUNN");
}

TEST(rna5, assign_char_member)
{
    rna5 t0;
    t0.assign_char('A');
    EXPECT_EQ(t0, rna5::A);
    t0.assign_char('C');
    EXPECT_EQ(t0, rna5::C);
    t0.assign_char('G');
    EXPECT_EQ(t0, rna5::G);
    t0.assign_char('T');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    t0.assign_char('U');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    t0.assign_char('N');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    t0.assign_char('a');
    EXPECT_EQ(t0, rna5::A);
    t0.assign_char('c');
    EXPECT_EQ(t0, rna5::C);
    t0.assign_char('g');
    EXPECT_EQ(t0, rna5::G);
    t0.assign_char('t');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    t0.assign_char('u');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    t0.assign_char('n');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    t0.assign_char('z');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);
    t0.assign_char('H');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);
    t0.assign_char('*');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    static_assert(std::is_same_v<decltype(t0.assign_char('C')), rna5 &>);
    EXPECT_EQ(t0.assign_char('C'), rna5::C);
}

TEST(rna5, assign_char_free)
{
    rna5 t0;
    assign_char(t0, 'A');
    EXPECT_EQ(t0, rna5::A);
    assign_char(t0, 'C');
    EXPECT_EQ(t0, rna5::C);
    assign_char(t0, 'G');
    EXPECT_EQ(t0, rna5::G);
    assign_char(t0, 'T');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    assign_char(t0, 'U');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    assign_char(t0, 'N');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    assign_char(t0, 'a');
    EXPECT_EQ(t0, rna5::A);
    assign_char(t0, 'c');
    EXPECT_EQ(t0, rna5::C);
    assign_char(t0, 'g');
    EXPECT_EQ(t0, rna5::G);
    assign_char(t0, 't');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    assign_char(t0, 'u');
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    assign_char(t0, 'n');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    assign_char(t0, 'z');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);
    assign_char(t0, 'H');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);
    assign_char(t0, '*');
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    static_assert(std::is_same_v<decltype(assign_char(t0, 'C')), rna5 &>);
    EXPECT_EQ(assign_char(t0, 'C'), rna5::C);
}

TEST(rna5, assign_rank_member)
{
    rna5 t0;
    t0.assign_rank(0);
    EXPECT_EQ(t0, rna5::A);
    t0.assign_rank(1);
    EXPECT_EQ(t0, rna5::C);
    t0.assign_rank(2);
    EXPECT_EQ(t0, rna5::G);
    t0.assign_rank(3);
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    t0.assign_rank(4);
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    static_assert(std::is_same_v<decltype(t0.assign_rank(2)), rna5 &>);
    EXPECT_EQ(t0.assign_rank(1), rna5::C);
}

TEST(rna5, assign_rank_free)
{
    rna5 t0;
    assign_rank(t0, 0);
    EXPECT_EQ(t0, rna5::A);
    assign_rank(t0, 1);
    EXPECT_EQ(t0, rna5::C);
    assign_rank(t0, 2);
    EXPECT_EQ(t0, rna5::G);
    assign_rank(t0, 3);
    EXPECT_EQ(t0, rna5::T);
    EXPECT_EQ(t0, rna5::U);
    assign_rank(t0, 4);
    EXPECT_EQ(t0, rna5::N);
    EXPECT_EQ(t0, rna5::UNKNOWN);

    static_assert(std::is_same_v<decltype(assign_rank(t0, 2)), rna5 &>);
    EXPECT_EQ(assign_rank(t0, 1), rna5::C);
}

// ------------------------------------------------------------------
// compatibility to dna5
// ------------------------------------------------------------------

TEST(rna5, to_dna5)
{
    rna5 t0{rna5::C};
    dna5 d0{t0};
    EXPECT_EQ(d0, dna5::C);

    d0 = rna5::U;
    EXPECT_EQ(d0, dna5::T);
}

TEST(rna5, from_dna5)
{
    dna5 d0{dna5::C};
    rna5 t0{d0};
    EXPECT_EQ(t0, rna5::C);

    t0 = dna5::T;
    EXPECT_EQ(t0, rna5::U);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(rna5_literals, vector)
{
    rna5_vector v;
    v.resize(5, rna5::A);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<rna5> w{rna5::A, rna5::C, rna5::G, rna5::T, rna5::U, rna5::N, rna5::UNKNOWN};
    EXPECT_EQ(w, "ACGUUNN"_rna5);
}

TEST(rna5_literals, basic_string)
{
    rna5_string v;
    v.resize(5, rna5::A);
    EXPECT_EQ(v, "AAAAA"_rna5s);

    std::basic_string<rna5, std::char_traits<rna5>> w{rna5::A, rna5::C, rna5::G, rna5::T, rna5::U, rna5::N,
                                                      rna5::UNKNOWN};
    EXPECT_EQ(w, "ACGUUNN"_rna5s);
}
