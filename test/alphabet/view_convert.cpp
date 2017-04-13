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

#include <sstream>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4_container.hpp>
#include <seqan3/alphabet/quality.hpp>
#include <seqan3/alphabet/view_convert.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(convert, to_char)
{
    dna4_vector vec = "ACTTTGATA"_dna4;

    auto v = vec | view::convert<char>;
    EXPECT_EQ(std::string{v}, "ACTTTGATA");

    std::vector<illumina18> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
    auto v2 = qvec | view::convert<char>;
    EXPECT_EQ(std::string{v2}, "!(&$(%?18");
}

TEST(convert, to_integral)
{
    dna4_vector vec = "ACTTTGATA"_dna4;

    std::vector<unsigned> v = vec | view::convert<unsigned>;
    std::vector<unsigned> const out{0,1,3,3,3,2,0,3,0};

    EXPECT_EQ(v, out);

    std::vector<illumina18> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
    std::vector<unsigned> v2 = qvec | view::convert<unsigned>;
    std::vector<unsigned> const out2{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
    EXPECT_EQ(v2, out2);
}

TEST(convert, default_impl)
{
    std::vector<dna4q> qcvec{{dna4::C, 0}, {dna4::A, 7}, {dna4::G, 5}, {dna4::T, 3}, {dna4::G, 7}, {dna4::A, 4},
                             {dna4::C, 30}, {dna4::T, 16}, {dna4::A, 23}};

    auto v6 = qcvec | view::convert<dna4> | view::convert<char>;
    EXPECT_EQ(std::string{v6}, "CAGTGACTA");
    auto v8 = qcvec | view::convert<illumina18> | view::convert<char>;
    EXPECT_EQ(std::string{v8}, "!(&$(%?18");
}
