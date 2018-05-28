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
#include <seqan3/alphabet/quality/phred68.hpp>

using namespace seqan3;

class phred68_fixture : public ::testing::Test
{
public:
    // the internal offset_phred
    const std::int8_t offset = -5;
};

// constructor
TEST(phred68_ctr, ctr)
{
    [[maybe_unused]] phred68 phred;
}

// default copy constructor
TEST(phred68_cp_ctr, cp_ctr)
{
    // in value_size range
    phred68 phred1{0};
    phred68 phred1_cp(phred1);
    phred68 phred2{67};
    phred68 phred2_cp(phred2);
}

// default destructor
TEST(phred68_des, des)
{
    phred68* phred_ptr = new phred68{};
    delete phred_ptr;
}

// cp by assignment
TEST(phred68_cp_ass, cp_ass)
{
    phred68 phred{0};
    [[maybe_unused]] phred68 phred2 = phred;
}

// internal ascii character offset (see https://en.wikipedia.org/wiki/FASTQ_format)
TEST(phred68_char_offset, const_offset)
{
    EXPECT_EQ(phred68::offset_char, ';');
}

// global and static quality alphabet size
TEST(phred68_alphabet_size, const_value_size)
{
    EXPECT_EQ(phred68::value_size, 68);
    EXPECT_EQ(alphabet_size_v<phred68>, 68);
}

// implicit value assignment
TEST(phred68_implicit_assign, implicit_assign)
{
    phred68 phred;
    phred = 19;
    // expect size unmodified
    EXPECT_EQ(phred68::value_size, 68);
    // newly assigned member, phred_score 19 has rank 19-(-5)
    EXPECT_EQ(phred._value, 24);
}

// test boolean ops
TEST(phred68_bool, compare)
{
    phred68 phred1{}, phred2{}, phred3{};
    phred1 = -3, phred2 = -3, phred3 = 0;
    EXPECT_TRUE(phred1 == phred2);
    EXPECT_TRUE(phred1 != phred3);
    EXPECT_TRUE(phred1 <= phred3);
    EXPECT_TRUE(phred1 < phred3);
    EXPECT_TRUE(phred3 >= phred1);
    EXPECT_TRUE(phred3 > phred1);
}

TEST_F(phred68_fixture, to_rank)
{
    phred68 phred;
    phred = 19 + offset;
    EXPECT_EQ(19, to_rank(phred));
    EXPECT_EQ(19, phred.to_rank());
    phred = 62;
    EXPECT_EQ(67, to_rank(phred));
    EXPECT_EQ(67, phred.to_rank());
}

// global assign_char operator
TEST(phred68_assign_char, assign_char)
{
    phred68 phred;
    phred = assign_char(phred, ';');
    EXPECT_EQ(0, to_rank(phred));
    phred = assign_char(phred, 'J');
    EXPECT_EQ(15, to_rank(phred));
}

// global and internal to_char
TEST_F(phred68_fixture, op_to_char)
{
    phred68 phred;
    phred = 2 + offset;
    EXPECT_EQ(to_char(phred), '=');
    EXPECT_EQ(phred.to_char(), '=');
    // change internal value
    phred = 62 + offset;
    EXPECT_EQ(to_char(phred), 'y');
    EXPECT_EQ(phred.to_char(), 'y');
}

TEST_F(phred68_fixture, assign_phred)
{
    phred68 phred;
    phred = 7 + offset;
    phred = assign_phred(phred, 9 + offset);
    [[maybe_unused]] seqan3::phred68::rank_type val = phred._value;
    EXPECT_EQ(9, to_rank(phred));
}

TEST(phred68_assign_rank, assign_rank)
{
    phred68 phred;
    phred = 7;
    phred = assign_rank(phred, 9);
    EXPECT_EQ(9, to_rank(phred));
}

TEST_F(phred68_fixture, to_phred)
{
    phred68 phred{}; // expect internal rank value 0 per default
    EXPECT_EQ(offset, to_phred(phred));
    phred = 39; // assign phred
    EXPECT_EQ(39, to_phred(phred));
    phred = 42;
    EXPECT_EQ(42, to_phred(phred));
}

TEST(phred68_cmp, cmp)
{
    // phred68 p{<num>} should never be called due to unimplemented range check
    phred68 phred1, phred2, phred3;
    // assign phred values [offset .. offset+68[
    phred1 = 7, phred2 = 11, phred3 = 62;

    EXPECT_LT(phred1, phred2);
    EXPECT_LE(phred1, phred2);
    EXPECT_LE(phred2, phred2);
    EXPECT_EQ(phred2, phred2);
    EXPECT_GE(phred2, phred2);
    EXPECT_GE(phred3, phred2);
    EXPECT_GT(phred3, phred2);
}
