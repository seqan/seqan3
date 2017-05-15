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
// Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/illumina18.hpp>

using namespace seqan3;

// constructor
TEST(illumina18_ctr, ctr)
{
    illumina18 illu;
}

// default copy constructor
TEST(illumina18_cp_ctr, cp_ctr)
{
    illumina18 illu;
    illumina18 illu2(illu);
}

// default destructor
TEST(illumina18_des, des)
{
    illumina18* illu_ptr = new illumina18{};
    delete illu_ptr;
}

// cp by assignment
TEST(illumina18_cp_ass, cp_ass)
{
    illumina18 illu;
    illumina18 illu2 = illu;
}

// phred score offset
TEST(illumina18_int_offset, int_offset)
{
    EXPECT_EQ(illumina18::offset_phred, 0);
}

// char offset
TEST(illumina18_char_offset, const_offset)
{
    EXPECT_EQ(illumina18::offset_char, '!');
}

// global and static quality alphabet size
TEST(illumina18_alphabet_size, const_value_size)
{
    EXPECT_EQ(illumina18::value_size, 42);
    EXPECT_EQ(alphabet_size_v<illumina18>, 42);
}

// implicit value assignment
TEST(illumina18_implicit_assign, implicit_assign)
{
    illumina18 illu;
    illu = 19;
    // expect size unmodified
    EXPECT_EQ(illumina18::value_size, 42);
    // newly assigned member
    EXPECT_EQ(illu.value, 19);
}

// char operator
TEST(illumina18_op_char, op_char)
{
    illumina18 illu;
    illu = 0;
    char c = char(illu);
    EXPECT_EQ(c, '!');
}

TEST(illumina18_assign_rank, assign_rank)
{
    illumina18 illu;
    illumina18 illu2 = assign_rank(illu, 1);
    EXPECT_EQ(1, to_rank(illu2));

    illu2 = illu.assign_rank(2);
    EXPECT_EQ(2, to_rank(illu2));
}

TEST(illumina18_to_rank, to_rank)
{
    illumina18 illu;
    illu = 19;
    EXPECT_EQ(19, to_rank(illu));
    EXPECT_EQ(19, illu.to_rank());
}

// global assign_char operator
TEST(illumina18_assign_char, assign_char)
{
    illumina18 illu;
    illu = assign_char(illu, '!');
    EXPECT_EQ(0, to_rank(illu));
}

// global and internal to_char
TEST(illumina18_op_tochar, op_to_char)
{
    illumina18 illu;
    illu = 2;
    EXPECT_EQ(to_char(illu), '#');
    EXPECT_EQ(illu.to_char(), '#');
    // change internal value
    illu = 41;
    EXPECT_EQ(to_char(illu), 'J');
    EXPECT_EQ(illu.to_char(), 'J');
}

TEST(illumina18_assign_phred, assign_phred)
{
    illumina18 illu;
    illu = 7;
    illu = assign_phred(illu, 9);
    seqan3::illumina18::rank_type val = illu.value;
    EXPECT_EQ(9, to_rank(illu));
}

TEST(illumina18_to_phred, to_phred)
{
    illumina18 illu{};
    EXPECT_EQ(0, to_phred(illu));
    illu = 39;
    EXPECT_EQ(39, to_rank(illu));
}

TEST(illumina18_cmp, cmp)
{
    illumina18 illu1{7};
    illumina18 illu2{11};
    illumina18 illu3{30};

    EXPECT_LT(illu1, illu2);
    EXPECT_LE(illu1, illu2);
    EXPECT_LE(illu2, illu2);
    EXPECT_EQ(illu2, illu2);
    EXPECT_GE(illu2, illu2);
    EXPECT_GE(illu3, illu2);
    EXPECT_GT(illu3, illu2);
}
