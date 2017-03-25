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

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gapped_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

TEST(gapped_alphabet_test, test_alphabet_concept)
{
    EXPECT_TRUE(alphabet_concept<gapped_alphabet<dna4>>);
}

TEST(gapped_alphabet_test, test_implicit_inner_type_compatibility)
{
    EXPECT_EQ(gapped_alphabet<dna4>{0}, gapped_alphabet<dna4>{}.from_integral(0));
}

TEST(gapped_alphabet_test, test_default_initialization)
{
    EXPECT_EQ(gapped_alphabet<dna4>{}, gapped_alphabet<dna4>{0});
}

TEST(gapped_alphabet_test, test_relations)
{
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('A'), gapped_alphabet<dna4>{}.from_char('A'));
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('a'), gapped_alphabet<dna4>{}.from_char('A'));
    EXPECT_NE(gapped_alphabet<dna4>{}.from_char('A'), gapped_alphabet<dna4>{}.from_char('C'));
    EXPECT_NE(gapped_alphabet<dna4>{}.from_char('A'), gapped_alphabet<dna4>{}.from_char('-'));
    EXPECT_LT(gapped_alphabet<dna4>{}.from_char('A'), gapped_alphabet<dna4>{}.from_char('C'));
    EXPECT_LE(gapped_alphabet<dna4>{}.from_char('C'), gapped_alphabet<dna4>{}.from_char('C'));
    EXPECT_LE(gapped_alphabet<dna4>{}.from_char('A'), gapped_alphabet<dna4>{}.from_char('C'));
    EXPECT_GT(gapped_alphabet<dna4>{}.from_char('T'), gapped_alphabet<dna4>{}.from_char('A'));
    EXPECT_GE(gapped_alphabet<dna4>{}.from_char('T'), gapped_alphabet<dna4>{}.from_char('T'));
    EXPECT_GE(gapped_alphabet<dna4>{}.from_char('T'), gapped_alphabet<dna4>{}.from_char('C'));
}

TEST(gapped_alphabet_test, test_stream_operator)
{
    std::stringstream ss;
    ss << gapped_alphabet<dna4>{0} << gapped_alphabet<dna4>{3} << gapped_alphabet<dna4>{2} << gapped_alphabet<dna4>{4} << gapped_alphabet<dna4>{1};
    EXPECT_EQ(ss.str(), "ATG-C");
}

TEST(gapped_alphabet_test, test_to_char)
{
    EXPECT_EQ(gapped_alphabet<dna4>{0}.to_char(), 'A');
    EXPECT_EQ(gapped_alphabet<dna4>{1}.to_char(), 'C');
    EXPECT_EQ(gapped_alphabet<dna4>{2}.to_char(), 'G');
    EXPECT_EQ(gapped_alphabet<dna4>{3}.to_char(), 'T');
    EXPECT_EQ(gapped_alphabet<dna4>{4}.to_char(), '-');
}

TEST( gapped_alphabet_test, test_from_char)
{
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('A'), gapped_alphabet<dna4>{0});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('C'), gapped_alphabet<dna4>{1});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('G'), gapped_alphabet<dna4>{2});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('T'), gapped_alphabet<dna4>{3});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('-'), gapped_alphabet<dna4>{4});
}

TEST(gapped_alphabet_test, test_to_integral)
{
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('A').to_integral(), 0);
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('C').to_integral(), 1);
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('G').to_integral(), 2);
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('T').to_integral(), 3);
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_char('-').to_integral(), 4);
}

TEST(gapped_alphabet_test, test_from_integral)
{
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_integral(0), gapped_alphabet<dna4>{0});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_integral(1), gapped_alphabet<dna4>{1});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_integral(2), gapped_alphabet<dna4>{2});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_integral(3), gapped_alphabet<dna4>{3});
    EXPECT_EQ(gapped_alphabet<dna4>{}.from_integral(4), gapped_alphabet<dna4>{4});
}
