// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>

using namespace seqan3;

// constructor
TEST(phred63_ctr, ctr)
{
    [[maybe_unused]] phred63 illu;
}

// default copy constructor
TEST(phred63_cp_ctr, cp_ctr)
{
    // in value_size range
    phred63 illu1{0};
    phred63 illu1_cp(illu1);
    phred63 illu2{52};
    phred63 illu2_cp(illu2);
}

// default destructor
TEST(phred63_des, des)
{
    phred63* illu_ptr = new phred63{};
    delete illu_ptr;
}

// cp by assignment
TEST(phred63_cp_ass, cp_ass)
{
    phred63 illu{0};
    [[maybe_unused]] phred63 illu2 = illu;
}

// char offset
TEST(phred63_char_offset, const_offset)
{
    EXPECT_EQ(phred63::offset_char, '!');
}

// global and static quality alphabet size
TEST(phred63_alphabet_size, const_value_size)
{
    EXPECT_EQ(phred63::value_size, 63);
    EXPECT_EQ(alphabet_size_v<phred63>, 63);
}

// implicit value assignment
TEST(phred63_implicit_assign, implicit_assign)
{
    phred63 illu;
    illu = 19;
    // expect size unmodified
    EXPECT_EQ(phred63::value_size, 63);
    // newly assigned member
    EXPECT_EQ(illu._value, 19);
}

TEST(phred63_to_rank, to_rank)
{
    phred63 illu;
    illu = 19;
    EXPECT_EQ(19, to_rank(illu));
    EXPECT_EQ(19, illu.to_rank());
    illu = 62;
    EXPECT_EQ(62, to_rank(illu));
    EXPECT_EQ(62, illu.to_rank());
}

// global assign_char operator
TEST(phred63_assign_char, assign_char)
{
    phred63 illu;
    illu = assign_char(illu, '!');
    EXPECT_EQ(0, to_rank(illu));
    illu = assign_char(illu, '1');
    EXPECT_EQ(16, to_rank(illu));
}

// global and internal to_char
TEST(phred63_op_tochar, op_to_char)
{
    phred63 illu;
    illu = 2;
    EXPECT_EQ(to_char(illu), '#');
    EXPECT_EQ(illu.to_char(), '#');
    // change internal value
    illu = 62;
    EXPECT_EQ(to_char(illu), '_');
    EXPECT_EQ(illu.to_char(), '_');
}

TEST(phred63_assign_phred, assign_phred)
{
    phred63 illu;
    illu = 7;
    illu = assign_phred(illu, 9);
    [[maybe_unused]] seqan3::phred63::rank_type val = illu._value;
    EXPECT_EQ(9, to_rank(illu));
}

TEST(phred63_to_phred, to_phred)
{
    phred63 illu;
    illu = 0;
    EXPECT_EQ(0, to_phred(illu));
    illu = 62;
    EXPECT_EQ(62, to_phred(illu));
}

TEST(phred63_cmp, cmp)
{
    // phred63 p{<num>} should never be called due to unimplemented range check
    phred63 illu1, illu2, illu3;
    illu1 = 7, illu2 = 11, illu3 = 59;

    EXPECT_LT(illu1, illu2);
    EXPECT_LE(illu1, illu2);
    EXPECT_LE(illu2, illu2);
    EXPECT_EQ(illu2, illu2);
    EXPECT_GE(illu2, illu2);
    EXPECT_GE(illu3, illu2);
    EXPECT_GT(illu3, illu2);
}
