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

#include <gtest/gtest.h>
#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

TEST(dna5_test, test_alphabet_concept)
{
    EXPECT_TRUE(alphabet_concept<dna5>);
}

TEST(dna5_test, test_default_initialization)
{
    EXPECT_EQ(dna5{dna5::A}, dna5{});
}

TEST(dna5_test, test_implicit_inner_type_compatibility)
{
    EXPECT_EQ(dna5{dna5::A}, dna5::A);
}

TEST(dna5_test, test_relations)
{
    EXPECT_LT(dna5{dna5::A}, dna5{dna5::C});
}

TEST(dna5_test, test_static_cast)
{
    EXPECT_EQ(static_cast<dna5::char_type>(dna5{dna5::C}), 'C');
}

TEST(dna5_test, test_stream_operator)
{
    std::stringstream ss;
    ss << dna5{dna5::A} << dna5{dna5::C} << dna5{dna5::G} << dna5{dna5::T} << dna5{dna5::N};
    EXPECT_EQ(ss.str(), "ACGTN");
}

TEST(dna5_test, test_to_char)
{
    EXPECT_EQ(dna5{dna5::A}.to_char(), 'A');
    EXPECT_EQ(dna5{dna5::C}.to_char(), 'C');
    EXPECT_EQ(dna5{dna5::G}.to_char(), 'G');
    EXPECT_EQ(dna5{dna5::T}.to_char(), 'T');
    EXPECT_EQ(dna5{dna5::N}.to_char(), 'N');
}

TEST( dna5_test, test_from_char)
{
    EXPECT_EQ(dna5{}.from_char('A'), dna5{dna5::A});
    EXPECT_EQ(dna5{}.from_char('C'), dna5{dna5::C});
    EXPECT_EQ(dna5{}.from_char('G'), dna5{dna5::G});
    EXPECT_EQ(dna5{}.from_char('T'), dna5{dna5::T});
    EXPECT_EQ(dna5{}.from_char('N'), dna5{dna5::N});

    EXPECT_EQ(dna5{}.from_char('a'), dna5{dna5::A});
    EXPECT_EQ(dna5{}.from_char('c'), dna5{dna5::C});
    EXPECT_EQ(dna5{}.from_char('g'), dna5{dna5::G});
    EXPECT_EQ(dna5{}.from_char('t'), dna5{dna5::T});
    EXPECT_EQ(dna5{}.from_char('n'), dna5{dna5::N});

    EXPECT_EQ(dna5{}.from_char('R'), dna5{dna5::N});
    EXPECT_EQ(dna5{}.from_char('x'), dna5{dna5::N});
    EXPECT_EQ(dna5{}.from_char('8'), dna5{dna5::N});
}

TEST(dna5_test, test_to_integral)
{
    EXPECT_EQ(dna5{dna5::A}.to_integral(), 0);
    EXPECT_EQ(dna5{dna5::C}.to_integral(), 1);
    EXPECT_EQ(dna5{dna5::G}.to_integral(), 2);
    EXPECT_EQ(dna5{dna5::T}.to_integral(), 3);
    EXPECT_EQ(dna5{dna5::N}.to_integral(), 4);
}

TEST(dna5_test, test_from_integral)
{
    EXPECT_EQ(dna5{}.from_integral(0), dna5{dna5::A});
    EXPECT_EQ(dna5{}.from_integral(1), dna5{dna5::C});
    EXPECT_EQ(dna5{}.from_integral(2), dna5{dna5::G});
    EXPECT_EQ(dna5{}.from_integral(3), dna5{dna5::T});
    EXPECT_EQ(dna5{}.from_integral(4), dna5{dna5::N});
}
