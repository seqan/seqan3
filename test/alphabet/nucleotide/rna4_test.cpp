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
// Author: Chenxu Pan <Chenxu Pan@fu-berlin.de>
// ==========================================================================
// Test cases for the biological rna4 alphabet.
// ==========================================================================

#include <seqan3/alphabet/nucleotide/rna4.hpp>
//#include "/home/cxpan/tmp/seqan3/seqan3/include/seqan3/alphabet/nucleotide/rna4.hpp"
#include <gtest/gtest.h>
//#include "/home/cxpan/tmp/seqan3/seqan3/vendor/googletest/googletest/include/gtest/gtest.h"
#include <sstream>

using namespace seqan3;

TEST(rna4_test, test_alphabet_concept)
{
    EXPECT_TRUE(alphabet_concept<rna4>);
}

TEST(rna4_test, test_default_initialization)
{
    EXPECT_EQ(rna4{rna4::A}, rna4{});
}

TEST(rna4_test, test_implicit_inner_type_compatibility)
{
    EXPECT_EQ(rna4{rna4::A}, rna4::A);
}

TEST(rna4_test, test_relations)
{
    EXPECT_LT(rna4{rna4::A}, rna4{rna4::C});
}

TEST(rna4_test, test_static_cast)
{
    EXPECT_EQ(static_cast<rna4::char_type>(rna4{rna4::C}), 'C');
}

TEST(rna4_test, test_stream_operator)
{
    std::stringstream ss;
    ss << rna4{rna4::A} << rna4{rna4::C} << rna4{rna4::G} << rna4{rna4::U};
    EXPECT_EQ(ss.str(), "ACGU");
}

TEST(rna4_test, test_to_char)
{
    EXPECT_EQ(rna4{rna4::A}.to_char(), 'A');
    EXPECT_EQ(rna4{rna4::C}.to_char(), 'C');
    EXPECT_EQ(rna4{rna4::G}.to_char(), 'G');
    EXPECT_EQ(rna4{rna4::T}.to_char(), 'U');
}

TEST( rna4_test, test_from_char)
{
    EXPECT_EQ(rna4{}.from_char('A'), rna4{rna4::A});
    EXPECT_EQ(rna4{}.from_char('C'), rna4{rna4::C});
    EXPECT_EQ(rna4{}.from_char('G'), rna4{rna4::G});
    EXPECT_EQ(rna4{}.from_char('U'), rna4{rna4::U});

    EXPECT_EQ(rna4{}.from_char('a'), rna4{rna4::A});
    EXPECT_EQ(rna4{}.from_char('c'), rna4{rna4::C});
    EXPECT_EQ(rna4{}.from_char('g'), rna4{rna4::G});
    EXPECT_EQ(rna4{}.from_char('u'), rna4{rna4::U});

    EXPECT_EQ(rna4{}.from_char('R'), rna4{rna4::A});
    EXPECT_EQ(rna4{}.from_char('x'), rna4{rna4::A});
    EXPECT_EQ(rna4{}.from_char('8'), rna4{rna4::A});
}

TEST(rna4_test, test_to_integral)
{
    EXPECT_EQ(rna4{rna4::A}.to_integral(), 0);
    EXPECT_EQ(rna4{rna4::C}.to_integral(), 1);
    EXPECT_EQ(rna4{rna4::G}.to_integral(), 2);
    EXPECT_EQ(rna4{rna4::U}.to_integral(), 3);
}

TEST(rna4_test, test_from_integral)
{
    EXPECT_EQ(rna4{}.from_integral(0), rna4{rna4::A});
    EXPECT_EQ(rna4{}.from_integral(1), rna4{rna4::C});
    EXPECT_EQ(rna4{}.from_integral(2), rna4{rna4::G});
    EXPECT_EQ(rna4{}.from_integral(3), rna4{rna4::U});
}
