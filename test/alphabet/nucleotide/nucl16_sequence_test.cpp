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
// Author: Sara Hetzel <sara.hetzel AT fu-berlin.de>
// ============================================================================

#include <array>

#include <gtest/gtest.h>
#include <seqan3/alphabet/nucleotide/nucl16.hpp>
#include <seqan3/alphabet/nucleotide/nucl16_sequence.hpp>

using namespace seqan3;

TEST(alphabet_nucleotides_nucl16_sequence_test, vector)
{
    nucl16_vector v{{nucl16::A}, {nucl16::C}, {nucl16::G}, {nucl16::T}};

    EXPECT_EQ(v[0], nucl16::A);
    EXPECT_EQ(v[1], nucl16::C);
    EXPECT_EQ(v[2], nucl16::G);
    EXPECT_EQ(v[3], nucl16::T);
}

TEST(alphabet_nucleotides_nucl16_sequence_test, string)
{
    nucl16_string s{{nucl16::A}, {nucl16::C}, {nucl16::G}, {nucl16::T}};
    EXPECT_EQ(s[0], nucl16::A);
    EXPECT_EQ(s[1], nucl16::C);
    EXPECT_EQ(s[2], nucl16::G);
    EXPECT_EQ(s[3], nucl16::T);
}
