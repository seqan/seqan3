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
// Author: David Heller <david.heller@fu-berlin.de>
// ==========================================================================
// Test cases for the biological dna5 alphabet.
// ==========================================================================

#include <seqan3/alphabet/nucleotide/dna5_sequence.hpp>
#include <gtest/gtest.h>
#include <sstream>
#include <vector>

using namespace seqan3;
using namespace seqan3::literal;

TEST(dna5_test, test_dna5_vector_operator)
{
    dna5_vector v;
    v.resize(5, dna5{dna5::A});
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<dna5> w {dna5{dna5::A}, dna5{dna5::C}, dna5{dna5::G}, dna5{dna5::T}, dna5{dna5::N}};
    EXPECT_EQ(w, "ACGTN"_dna5);
}

TEST(dna5_test, test_dna5_string_operator)
{
    dna5_string v;
    v.resize(5, dna5{dna5::A});
    EXPECT_EQ(v, "AAAAA"_dna5s);

    std::basic_string<dna5, std::char_traits<dna5>> w {dna5{dna5::A}, dna5{dna5::C}, dna5{dna5::G},
                                                       dna5{dna5::T}, dna5{dna5::N}};
    EXPECT_EQ(w, "ACGTN"_dna5s);
}
