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
// Author: Chenxu Pan <chenxu.pan@fu-berlin.de>
// ==========================================================================
// Test cases for the biological rna4 alphabet.
// ==========================================================================

#include <seqan3/alphabet/nucleotide/rna4_sequence.hpp>
#include <gtest/gtest.h>
#include <sstream>
#include <vector>

using namespace seqan3;
using namespace seqan3::literal;

TEST(rna4_test, test_rna4_vector_operator)
{
    rna4_vector v;
    v.resize(4, rna4{rna4::A});
    EXPECT_EQ(v, "AAAA"_rna4);

    std::vector<rna4> w {rna4{rna4::A}, rna4{rna4::C}, rna4{rna4::G}, rna4{rna4::U}};
    EXPECT_EQ(w, "ACGU"_rna4);
}

TEST(rna4_test, test_rna4_string_operator)
{
    rna4_string v;
    v.resize(4, rna4{rna4::A});
    EXPECT_EQ(v, "AAAA"_rna4);

    std::basic_string<rna4, std::char_traits<rna4>> w {rna4{rna4::A}, rna4{rna4::C}, rna4{rna4::G},
                                                       rna4{rna4::U}};
    EXPECT_EQ(w, "ACGU"_rna4);
}
