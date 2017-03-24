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
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// ==========================================================================

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/union_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

TEST(union_alphabet_test, default_ctr)
{
    using alphabet_t = union_alphabet<dna4, gap>;
    alphabet_t{};
}

// TODO: othe constructurs

TEST(union_alphabet_test, assignment_union_base_types)
{
    // TODO:
    // using alphabet_t = union_alphabet<dna4, dna5, gap>;
    // alphabet_t letter1 = dna4{};
    // alphabet_t letter2 = dna5{};
    // alphabet_t letter3 = gap{};
}

TEST(union_alphabet_test, single_union)
{
    using alphabet_t = union_alphabet<dna4>;
    alphabet_t{};
}

TEST(union_alphabet_test, integral_type)
{
    using alphabet1_t = union_alphabet<dna4, dna5, gap>;
    using alphabet2_t = union_alphabet<gap, dna5, dna4>;
    using alphabet3_t = union_alphabet<gap>;

    constexpr auto expect1 = std::is_same_v<alphabet1_t::integral_type, uint8_t>;
    EXPECT_TRUE(expect1);

    constexpr auto expect2 = std::is_same_v<alphabet2_t::integral_type, uint8_t>;
    EXPECT_TRUE(expect2);

    constexpr auto expect3 = std::is_same_v<alphabet3_t::integral_type, bool>;
    EXPECT_TRUE(expect3);
}

// todo: constructors

TEST(union_alphabet_test, from_and_to_integral)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    alphabet_t letter{};

    EXPECT_EQ(letter.from_integral(0).to_integral(), 0);
    EXPECT_EQ(letter.from_integral(1).to_integral(), 1);
    EXPECT_EQ(letter.from_integral(2).to_integral(), 2);
    EXPECT_EQ(letter.from_integral(3).to_integral(), 3);
    EXPECT_EQ(letter.from_integral(4).to_integral(), 4);
    EXPECT_EQ(letter.from_integral(5).to_integral(), 5);
    EXPECT_EQ(letter.from_integral(6).to_integral(), 6);
    EXPECT_EQ(letter.from_integral(7).to_integral(), 7);
    EXPECT_EQ(letter.from_integral(8).to_integral(), 8);
    EXPECT_EQ(letter.from_integral(9).to_integral(), 9);
}

TEST(union_alphabet_test, to_char)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    alphabet_t letter{};

    EXPECT_EQ(letter.from_integral(0).to_char(), 'A');
    EXPECT_EQ(letter.from_integral(1).to_char(), 'C');
    EXPECT_EQ(letter.from_integral(2).to_char(), 'G');
    EXPECT_EQ(letter.from_integral(3).to_char(), 'T');
    EXPECT_EQ(letter.from_integral(4).to_char(), 'A');
    EXPECT_EQ(letter.from_integral(5).to_char(), 'C');
    EXPECT_EQ(letter.from_integral(6).to_char(), 'G');
    EXPECT_EQ(letter.from_integral(7).to_char(), 'T');
    EXPECT_EQ(letter.from_integral(8).to_char(), 'N');
    EXPECT_EQ(letter.from_integral(9).to_char(), '-');

    EXPECT_EQ(letter.from_integral(10).to_char(), static_cast<char>(0));
}

TEST(union_alphabet_test, relations)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    alphabet_t letter1{};
    alphabet_t letter2{};

    EXPECT_EQ(letter1.from_integral(0), letter2.from_integral(0));
    EXPECT_EQ(letter1.from_integral(5), letter2.from_integral(5));
    EXPECT_NE(letter1.from_integral(1), letter2.from_integral(5));
    EXPECT_GT(letter1.from_integral(2), letter2.from_integral(1));
    EXPECT_GE(letter1.from_integral(2), letter2.from_integral(1));
    EXPECT_GE(letter1.from_integral(2), letter2.from_integral(2));
    EXPECT_LT(letter1.from_integral(1), letter2.from_integral(2));
    EXPECT_LE(letter1.from_integral(1), letter2.from_integral(2));
    EXPECT_LE(letter1.from_integral(2), letter2.from_integral(2));
}
