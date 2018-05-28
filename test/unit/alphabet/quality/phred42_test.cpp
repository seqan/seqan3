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
#include <seqan3/alphabet/quality/phred42.hpp>

using namespace seqan3;

// constructor
TEST(phred42_ctr, ctr)
{
    [[maybe_unused]] phred42 phred;
}

// default copy constructor
TEST(phred42_cp_ctr, cp_ctr)
{
    // in value_size range
    phred42 phred1;
    phred1 = 3;
    phred42 phred1_cp(phred1);
    // in max_value_size range
    phred42 phred2;
    phred2 = 52;
    phred42 phred2_cp(phred2);
}

// default destructor
TEST(phred42_des, des)
{
    phred42* phred_ptr = new phred42{};
    delete phred_ptr;
}

// cp by assignment
TEST(phred42_cp_ass, cp_ass)
{
    phred42 phred;
    phred = 0;
    [[maybe_unused]] phred42 phred2 = phred;
}

// char offset
TEST(phred42_char_offset, const_offset)
{
    EXPECT_EQ(phred42::offset_char, '!');
}

// global and static quality alphabet size
TEST(phred42_alphabet_size, const_value_size)
{
    EXPECT_EQ(phred42::value_size, 42);
    EXPECT_EQ(alphabet_size_v<phred42>, 42);
}

// implicit value assignment
TEST(phred42_implicit_assign, implicit_assign)
{
    phred42 phred;
    phred = 19;
    // expect size unmodified
    EXPECT_EQ(phred42::value_size, 42);
    // newly assigned member
    EXPECT_EQ(phred._value, 19);
    // value in [value_size .. max_value_size[
    phred42 phred2;
    phred2 = 49;
    EXPECT_EQ(phred2._value, 41);
}

TEST(phred42_to_rank, to_rank)
{
    phred42 phred;
    phred = 19;
    EXPECT_EQ(19, to_rank(phred));
    EXPECT_EQ(19, phred.to_rank());
    // in [value_size .. max_value_size[
    phred = 47;
    EXPECT_EQ(41, to_rank(phred));
    EXPECT_EQ(41, phred.to_rank());
}

// global assign_char operator
TEST(phred42_assign_char, assign_char)
{
    phred42 phred;
    phred = assign_char(phred, '!');
    EXPECT_EQ(0, to_rank(phred));
    // in [value_size .. max_value_size[
    phred = assign_char(phred, 'O');
    EXPECT_EQ(41, to_rank(phred));
}

// global and internal to_char
TEST(phred42_op_tochar, op_to_char)
{
    phred42 phred;
    phred = 2;
    EXPECT_EQ(to_char(phred), '#');
    EXPECT_EQ(phred.to_char(), '#');
    // change internal value
    phred = 41;
    EXPECT_EQ(to_char(phred), 'J');
    EXPECT_EQ(phred.to_char(), 'J');
    phred = 42;
    EXPECT_EQ(to_char(phred), 'J');
    EXPECT_EQ(phred.to_char(), 'J');
}

TEST(phred42_assign_phred, assign_phred)
{
    phred42 phred;
    phred = 7;
    phred = assign_phred(phred, 9);
    [[maybe_unused]] seqan3::phred42::rank_type val = phred._value;
    EXPECT_EQ(9, to_rank(phred));
    phred = assign_phred(phred, 43);
    EXPECT_EQ(41, to_rank(phred));
}

TEST(phred42_to_phred, to_phred)
{
    phred42 phred;
    phred = 0;
    EXPECT_EQ(0, to_phred(phred));
    phred = 39;
    EXPECT_EQ(39, to_rank(phred));
    phred = 41;
    EXPECT_EQ(41, to_rank(phred));
    // in [value_size .. max_value_size[
    phred = 42;
    EXPECT_EQ(41, to_rank(phred));
}

TEST(phred42_cmp, cmp)
{
    // phred42 p{<num>} should never be called due to unimplemented range check
    phred42 phred1, phred2, phred3, phred4, phred5;
    phred1 = 7, phred2 = 11, phred3 = 30, phred4 = 41, phred5 = 43;
    EXPECT_LT(phred1, phred2);
    EXPECT_LE(phred1, phred2);
    EXPECT_LE(phred2, phred2);
    EXPECT_EQ(phred2, phred2);
    EXPECT_GE(phred2, phred2);
    EXPECT_GE(phred3, phred2);
    EXPECT_GT(phred3, phred2);
    EXPECT_EQ(phred4, phred5);
}
