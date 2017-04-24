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

#include <seqan3/alphabet/nucleotide/nucl16.hpp>

using namespace seqan3;
using namespace seqan3::literal;

constexpr char capital_char[] = {'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y'};
constexpr char lower_char[] = {'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 't', 'u', 'v', 'w', 'y'};
constexpr nucl16 all_nucl16[] = {nucl16::A, nucl16::B, nucl16::C, nucl16::D, nucl16::G, nucl16::H, nucl16::K,
                                 nucl16::M, nucl16::N, nucl16::R, nucl16::S, nucl16::T, nucl16::U, nucl16::V,
                                 nucl16::W, nucl16::Y};

// default/zero construction
TEST(nucl16, ctr)
{
    nucl16 t1;
}

// zero initialization
TEST(nucl16, zro)
{
    nucl16 t0{};
    EXPECT_EQ(t0, nucl16::A);
}

// copy construction
TEST(nucl16, cp_ctr)
{
    nucl16 t1{nucl16::C};
    nucl16 t2{t1};
    nucl16 t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(nucl16, mv_ctr)
{
    nucl16 t0{nucl16::C};
    nucl16 t1{nucl16::C};
    nucl16 t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    nucl16 t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(nucl16, cp_assgn)
{
    nucl16 t1{nucl16::C};
    nucl16 t2;
    nucl16 t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(nucl16, mv_assgn)
{
    nucl16 t0{nucl16::C};
    nucl16 t1{nucl16::C};
    nucl16 t2;
    nucl16 t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(nucl16, swap)
{
    nucl16 t0{nucl16::C};
    nucl16 t1{nucl16::C};
    nucl16 t2{};
    nucl16 t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// comparison
TEST(nucl16, cmp)
{
    nucl16 t0{nucl16::A};
    nucl16 t1{nucl16::C};
    nucl16 t2{nucl16::G};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

TEST(nucl16, to_char_member)
{
    for (size_t i = 0; i < 16; ++i)
        EXPECT_EQ(all_nucl16[i].to_char(), capital_char[i]);
}

TEST(nucl16, to_char_free)
{
    for (size_t i = 0; i < 16; ++i)
        EXPECT_EQ(to_char(all_nucl16[i]), capital_char[i]);
}

TEST(nucl16, to_rank_member)
{
    for (size_t i = 0; i < 16; ++i)
        EXPECT_EQ(all_nucl16[i].to_rank(), i);
}

TEST(nucl16, to_rank_free)
{
    for (size_t i = 0; i < 16; ++i)
        EXPECT_EQ(all_nucl16[i].to_rank(), i);
}

TEST(nucl16, stream_operator)
{
    std::stringstream ss;
    ss << nucl16::A << nucl16::C << nucl16::G << nucl16::T << nucl16::U << nucl16::N << nucl16::UNKNOWN;
    EXPECT_EQ(ss.str(), "ACGTUNN");
}

TEST(nucl16, assign_char_member)
{
    nucl16 t0;
    for (size_t i = 0; i < 16; ++i)
    {
        t0.assign_char(lower_char[i]);
        EXPECT_EQ(t0, all_nucl16[i]);
        t0.assign_char(capital_char[i]);
        EXPECT_EQ(t0, all_nucl16[i]);
    }

    t0.assign_char('z');
    EXPECT_EQ(t0, nucl16::N);
    EXPECT_EQ(t0, nucl16::UNKNOWN);
    t0.assign_char('*');
    EXPECT_EQ(t0, nucl16::N);
    EXPECT_EQ(t0, nucl16::UNKNOWN);

    static_assert(std::is_same_v<decltype(t0.assign_char('C')), nucl16 &>);
    EXPECT_EQ(t0.assign_char('C'), nucl16::C);
}

TEST(nucl16, assign_char_free)
{
    nucl16 t0;
    for (size_t i = 0; i < 16; ++i)
    {
        assign_char(t0, lower_char[i]);
        EXPECT_EQ(t0, all_nucl16[i]);
        assign_char(t0, capital_char[i]);
        EXPECT_EQ(t0, all_nucl16[i]);
    }

    assign_char(t0, 'z');
    EXPECT_EQ(t0, nucl16::N);
    EXPECT_EQ(t0, nucl16::UNKNOWN);
    assign_char(t0, '*');
    EXPECT_EQ(t0, nucl16::N);
    EXPECT_EQ(t0, nucl16::UNKNOWN);

    static_assert(std::is_same_v<decltype(assign_char(t0, 'C')), nucl16 &>);
    EXPECT_EQ(assign_char(t0, 'C'), nucl16::C);
}

TEST(nucl16, assign_rank_member)
{
    nucl16 t0;
    for (size_t i = 0; i < 16; ++i)
    {
        t0.assign_rank(i);
        EXPECT_EQ(t0, all_nucl16[i]);
    }

    static_assert(std::is_same_v<decltype(t0.assign_rank(2)), nucl16 &>);
    EXPECT_EQ(t0.assign_rank(2), nucl16::C);
}

TEST(nucl16, assign_rank_free)
{
    nucl16 t0;
    for (size_t i = 0; i < 16; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(t0, all_nucl16[i]);
    }

    static_assert(std::is_same_v<decltype(assign_rank(t0, 2)), nucl16 &>);
    EXPECT_EQ(assign_rank(t0, 2), nucl16::C);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(nucl16_literals, vector)
{
    nucl16_vector v;
    v.resize(5, nucl16::A);
    EXPECT_EQ(v, "AAAAA"_nucl16);

    std::vector<nucl16> w{nucl16::A, nucl16::C, nucl16::G, nucl16::T, nucl16::U, nucl16::N, nucl16::UNKNOWN};
    EXPECT_EQ(w, "ACGTUNN"_nucl16);
}

TEST(nucl16_literals, basic_string)
{
    nucl16_string v;
    v.resize(5, nucl16::A);
    EXPECT_EQ(v, "AAAAA"_nucl16s);

    std::basic_string<nucl16, std::char_traits<nucl16>> w{nucl16::A, nucl16::C, nucl16::G, nucl16::T, nucl16::U, nucl16::N,
                                                      nucl16::UNKNOWN};
    EXPECT_EQ(w, "ACGTUNN"_nucl16s);
}
