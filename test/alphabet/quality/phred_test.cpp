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

#include "../../../include/seqan3/alphabet/quality.hpp"
#include "../../../include/seqan3/alphabet/quality/phred.hpp"

using namespace seqan3;

TEST(illumina18_constructor, default_ctr)
{
    illumina18 illu;
    EXPECT_EQ(illumina18::value_size, 42);
    EXPECT_EQ(illu.offset_char, '!');
}

TEST(illumina18_implicit_assign, implicit_assign)
{
    illumina18 illu;
    illu = 19;
    // expect size unmodified
    EXPECT_EQ(illumina18::value_size, 42);
    // newly assigned member
    EXPECT_EQ(illu.value, 19);
}

TEST(illumina18_op_char, op_char)
{
    illumina18 illu;
    illu = 0;
    char c = char(illu);
    EXPECT_EQ(c, '!');
}

TEST(illumina18_op_tochar, op_tochar)
{
    illumina18 illu;
    illu = 0;
    char c = illu.to_char();
    EXPECT_EQ(c, '!');
    illu = 41;
    c = illu.to_char();
    EXPECT_EQ(c, 'J');
}

TEST(illumina18_from_char, from_char)
{
    illumina18 illu;
    illu = illu.from_char('!');
    EXPECT_EQ(0, illu.value);
}

TEST(illumina18_to_integral, to_integral)
{
    illumina18 illu;
    illu = 0;
    illumina18 illu2 = illu.from_char('!');
    EXPECT_EQ(0, illu2.value);
}

TEST(illumina18_from_integral, from_integral)
{
    illumina18 illu;
    illu = 0;
    illumina18 illu2 = illu.from_integral(1);
    EXPECT_EQ(1, illu2.value);
}

TEST(illumina18_from_phred, from_phred)
{
    illumina18 illu;
    illumina18 illu2 = illu.from_phred(5);
    EXPECT_EQ(5, illu2.value);
}

TEST(illumina18_to_phred, from_phred)
{
    illumina18 illu;
    illumina18 illu2 = illu.from_phred(5);
    EXPECT_EQ(5, illu2.value);
}
