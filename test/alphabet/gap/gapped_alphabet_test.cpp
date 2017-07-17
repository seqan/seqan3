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

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gapped_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

TEST(gapped_alphabet_test, default_constructor)
{
    using alphabet_t = gapped_alphabet<dna4>;

    constexpr auto letter1 = alphabet_t{};
    EXPECT_EQ(letter1._value, 0);
}

TEST(gapped_alphabet_test, initialize_from_component_alphabet)
{
    using alphabet_t = gapped_alphabet<dna4>;

    constexpr alphabet_t letter0{dna4::A};
    constexpr alphabet_t letter1 = dna4::C;
    constexpr alphabet_t letter2 = {dna4::G};
    constexpr alphabet_t letter3 = static_cast<alphabet_t>(dna4::T);

    alphabet_t letter4{dna4::A};
    alphabet_t letter5 = dna4::C;
    alphabet_t letter6 = {dna4::G};
    alphabet_t letter7 = static_cast<alphabet_t>(dna4::T);

    constexpr alphabet_t letter8{gap::GAP}; // letter3 = dna4::T; does not work
    alphabet_t letter9{gap::GAP};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 0);
    EXPECT_EQ(letter5.to_rank(), 1);
    EXPECT_EQ(letter6.to_rank(), 2);
    EXPECT_EQ(letter7.to_rank(), 3);
    EXPECT_EQ(letter8.to_rank(), 4);
    EXPECT_EQ(letter9.to_rank(), 4);
}

TEST(gapped_alphabet_test, assign_from_component_alphabet)
{
    using alphabet_t = gapped_alphabet<dna4>;
    alphabet_t letter{};

    letter = dna4::A;
    EXPECT_EQ(letter.to_rank(), 0);

    letter = {dna4::C}; // letter = {dna4::C}; does not work
    EXPECT_EQ(letter.to_rank(), 1);

    letter = static_cast<alphabet_t>(dna4::G);
    EXPECT_EQ(letter.to_rank(), 2);

    letter = {static_cast<alphabet_t>(dna4::T)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = gap::GAP;
    EXPECT_EQ(letter.to_rank(), 4);
}

TEST(gapped_alphabet_test, copy_constructor)
{
    using alphabet_t = gapped_alphabet<dna4>;
    constexpr alphabet_t letter1{dna4::T};
    alphabet_t letter2{letter1};

    EXPECT_EQ(letter1._value, 3);
    EXPECT_EQ(letter2._value, 3);
}

TEST(gapped_alphabet_test, move_constructor)
{
    using alphabet_t = gapped_alphabet<dna4>;
    alphabet_t letter1{dna4::G};
    alphabet_t letter2{std::move(letter1)};

    EXPECT_EQ(letter2._value, 2);
}

TEST(gapped_alphabet_test, copy_assignment)
{
    using alphabet_t = gapped_alphabet<dna4>;
    constexpr alphabet_t letter1{dna4::T};
    alphabet_t letter2, letter3;
    letter2 = letter1;
    letter3 = {letter1};

    EXPECT_EQ(letter1.to_rank(), 3);
    EXPECT_EQ(letter2.to_rank(), 3);
    EXPECT_EQ(letter3.to_rank(), 3);
}

TEST(gapped_alphabet_test, move_assignment)
{
    using alphabet_t = gapped_alphabet<dna4>;
    alphabet_t letter1 = dna4::G;
    alphabet_t letter2{std::move(letter1)};

    EXPECT_EQ(letter2.to_rank(), 2);
}

TEST(gapped_alphabet_test, fulfills_concepts)
{
    using alphabet_t = gapped_alphabet<dna4>;
    static_assert(std::is_pod_v<alphabet_t>);
    EXPECT_TRUE(alphabet_concept<alphabet_t>);
}

TEST( gapped_alphabet_test, assign_char)
{
    using alphabet_t = gapped_alphabet<dna4>;

    auto letter = alphabet_t{};
    auto letterA = alphabet_t{dna4::A};
    auto letterC = alphabet_t{dna4::C};
    auto letterG = alphabet_t{dna4::G};
    auto letterT = alphabet_t{dna4::T};
    auto letterD = alphabet_t{gap::GAP};

    EXPECT_EQ(letter.assign_char('A'), letterA);
    EXPECT_EQ(letter.assign_char('C'), letterC);
    EXPECT_EQ(letter.assign_char('G'), letterG);
    EXPECT_EQ(letter.assign_char('T'), letterT);
    EXPECT_EQ(letter.assign_char('-'), letterD);
}

TEST(gapped_alphabet_test, to_char)
{
    using alphabet_t = gapped_alphabet<dna4>;

    auto letterA = alphabet_t{dna4::A};
    auto letterC = alphabet_t{dna4::C};
    auto letterG = alphabet_t{dna4::G};
    auto letterT = alphabet_t{dna4::T};
    auto letterD = alphabet_t{gap::GAP};

    EXPECT_EQ(letterA.to_char(), 'A');
    EXPECT_EQ(letterC.to_char(), 'C');
    EXPECT_EQ(letterG.to_char(), 'G');
    EXPECT_EQ(letterT.to_char(), 'T');
    EXPECT_EQ(letterD.to_char(), '-');
}

TEST(gapped_alphabet_test, to_rank)
{
    using alphabet_t = gapped_alphabet<dna4>;
    auto letter = alphabet_t{};

    EXPECT_EQ(letter.assign_char('A').to_rank(), 0);
    EXPECT_EQ(letter.assign_char('C').to_rank(), 1);
    EXPECT_EQ(letter.assign_char('G').to_rank(), 2);
    EXPECT_EQ(letter.assign_char('T').to_rank(), 3);
    EXPECT_EQ(letter.assign_char('-').to_rank(), 4);
}

TEST(gapped_alphabet_test, assign_rank)
{
    using alphabet_t = gapped_alphabet<dna4>;

    auto letter = alphabet_t{};
    auto letterA = alphabet_t{dna4::A};
    auto letterC = alphabet_t{dna4::C};
    auto letterG = alphabet_t{dna4::G};
    auto letterT = alphabet_t{dna4::T};
    auto letterD = alphabet_t{gap::GAP};

    EXPECT_EQ(letter.assign_rank(0), letterA);
    EXPECT_EQ(letter.assign_rank(1), letterC);
    EXPECT_EQ(letter.assign_rank(2), letterG);
    EXPECT_EQ(letter.assign_rank(3), letterT);
    EXPECT_EQ(letter.assign_rank(4), letterD);
}

TEST(gapped_alphabet_test, relations)
{
    using alphabet_t = gapped_alphabet<dna4>;
    auto letter1 = alphabet_t{};
    auto letter2 = alphabet_t{};

    EXPECT_EQ(letter1.assign_char('A'), letter2.assign_char('A'));
    EXPECT_EQ(letter1.assign_char('a'), letter2.assign_char('A'));
    EXPECT_NE(letter1.assign_char('A'), letter2.assign_char('C'));
    EXPECT_NE(letter1.assign_char('A'), letter2.assign_char('-'));
    EXPECT_LT(letter1.assign_char('A'), letter2.assign_char('C'));
    EXPECT_LE(letter1.assign_char('C'), letter2.assign_char('C'));
    EXPECT_LE(letter1.assign_char('A'), letter2.assign_char('C'));
    EXPECT_GT(letter1.assign_char('T'), letter2.assign_char('A'));
    EXPECT_GE(letter1.assign_char('T'), letter2.assign_char('T'));
    EXPECT_GE(letter1.assign_char('T'), letter2.assign_char('C'));
}

TEST(gapped_alphabet_test, stream_operator)
{
    using alphabet_t = gapped_alphabet<dna4>;

    auto letterA = alphabet_t{dna4::A};
    auto letterC = alphabet_t{dna4::C};
    auto letterG = alphabet_t{dna4::G};
    auto letterT = alphabet_t{dna4::T};
    auto letterD = alphabet_t{gap::GAP};

    std::stringstream ss;
    ss << letterA << letterT << letterG << letterD << letterC;
    EXPECT_EQ(ss.str(), "ATG-C");
}
