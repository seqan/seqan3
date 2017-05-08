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

#include <gtest/gtest.h>
#include <sstream>

#include <seqan3/alphabet/nucleotide.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class nucleotide : public ::testing::Test
{};

// add all alphabets from the nucleotide sub module here
using nucleotide_types = ::testing::Types<dna4, dna5, rna4, rna5, nucl16>;

TYPED_TEST_CASE(nucleotide, nucleotide_types);

TYPED_TEST(nucleotide, static_members)
{
    TypeParam::A;
    TypeParam::C;
    TypeParam::G;
    TypeParam::T;
    TypeParam::U;
    TypeParam::UNKNOWN;
}

TYPED_TEST(nucleotide, assign_char)
{
    EXPECT_EQ((assign_char(TypeParam{}, 'A')), TypeParam::A);
    EXPECT_EQ((assign_char(TypeParam{}, 'C')), TypeParam::C);
    EXPECT_EQ((assign_char(TypeParam{}, 'G')), TypeParam::G);
    EXPECT_EQ((assign_char(TypeParam{}, 'T')), TypeParam::T);
    EXPECT_EQ((assign_char(TypeParam{}, 'U')), TypeParam::U);
    EXPECT_EQ((assign_char(TypeParam{}, '!')), TypeParam::UNKNOWN);
}

TYPED_TEST(nucleotide, to_char)
{
    EXPECT_EQ(to_char(TypeParam::A), 'A');
    EXPECT_EQ(to_char(TypeParam::C), 'C');
    EXPECT_EQ(to_char(TypeParam::G), 'G');

    if constexpr (std::is_same_v<TypeParam, rna4> || std::is_same_v<TypeParam, rna5>)
    {
        EXPECT_EQ(to_char(TypeParam::U), 'U');
        EXPECT_EQ(to_char(TypeParam::T), 'U');
    }
    else if constexpr (std::is_same_v<TypeParam, dna4> || std::is_same_v<TypeParam, dna5>)
    {
        EXPECT_EQ(to_char(TypeParam::U), 'T');
        EXPECT_EQ(to_char(TypeParam::T), 'T');
    } else
    {
        EXPECT_EQ(to_char(TypeParam::U), 'U');
        EXPECT_EQ(to_char(TypeParam::T), 'T');
    }

    if constexpr (std::is_same_v<TypeParam, dna4> || std::is_same_v<TypeParam, rna4>)
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'A');
    else
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'N');
}

TYPED_TEST(nucleotide, stream_operator)
{
    std::stringstream ss;
    ss << TypeParam::A << TypeParam::C << TypeParam::G;
    EXPECT_EQ(ss.str(), "ACG");
}

TYPED_TEST(nucleotide, concept)
{
    EXPECT_TRUE(nucleotide_concept<TypeParam>);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(dna4_literals, vector)
{
    dna4_vector v;
    v.resize(5, dna4::A);
    EXPECT_EQ(v, "AAAAA"_dna4);

    std::vector<dna4> w{dna4::A, dna4::C, dna4::G, dna4::T, dna4::U, dna4::UNKNOWN};
    EXPECT_EQ(w, "ACGTTA"_dna4);
}

TEST(dna4_literals, basic_string)
{
    dna4_string v;
    v.resize(5, dna4::A);
    EXPECT_EQ(v, "AAAAA"_dna4s);

    std::basic_string<dna4, std::char_traits<dna4>> w{dna4::A, dna4::C, dna4::G, dna4::T, dna4::U, dna4::UNKNOWN};
    EXPECT_EQ(w, "ACGTTA"_dna4s);
}

TEST(dna5_literals, vector)
{
    dna5_vector v;
    v.resize(5, dna5::A);
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<dna5> w{dna5::A, dna5::C, dna5::G, dna5::T, dna5::U, dna5::N, dna5::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_dna5);
}

TEST(dna5_literals, basic_string)
{
    dna5_string v;
    v.resize(5, dna5::A);
    EXPECT_EQ(v, "AAAAA"_dna5s);

    std::basic_string<dna5, std::char_traits<dna5>> w{dna5::A, dna5::C, dna5::G, dna5::T, dna5::U, dna5::N,
                                                      dna5::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_dna5s);
}

TEST(rna4_literals, vector)
{
    rna4_vector v;
    v.resize(5, rna4::A);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<rna4> w{rna4::A, rna4::C, rna4::G, rna4::T, rna4::U, rna4::UNKNOWN};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

TEST(rna4_literals, basic_string)
{
    rna4_string v;
    v.resize(5, rna4::A);
    EXPECT_EQ(v, "AAAAA"_rna4s);

    std::basic_string<rna4, std::char_traits<rna4>> w{rna4::A, rna4::C, rna4::G, rna4::T, rna4::U, rna4::UNKNOWN};
    EXPECT_EQ(w, "ACGUUA"_rna4s);
}

TEST(rna5_literals, vector)
{
    rna5_vector v;
    v.resize(5, rna5::A);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<rna5> w{rna5::A, rna5::C, rna5::G, rna5::T, rna5::U, rna5::N, rna5::UNKNOWN};
    EXPECT_EQ(w, "ACGUUNN"_rna5);
}

TEST(rna5_literals, basic_string)
{
    rna5_string v;
    v.resize(5, rna5::A);
    EXPECT_EQ(v, "AAAAA"_rna5s);

    std::basic_string<rna5, std::char_traits<rna5>> w{rna5::A, rna5::C, rna5::G, rna5::T, rna5::U, rna5::N,
                                                      rna5::UNKNOWN};
    EXPECT_EQ(w, "ACGUUNN"_rna5s);
}

TEST(nucl16_literals, vector)
{
    nucl16_vector v;
    v.resize(5, nucl16::A);
    EXPECT_EQ(v, "AAAAA"_nucl16);

    std::vector<nucl16> w{nucl16::A, nucl16::C, nucl16::G, nucl16::T, nucl16::U, nucl16::N, nucl16::UNKNOWN};
    EXPECT_EQ(w, "ACGTUNN"_nucl16);
}

TEST(nucl16_literals, basic_string)
{
    nucl16_string v;
    v.resize(5, nucl16::A);
    EXPECT_EQ(v, "AAAAA"_nucl16s);

    std::basic_string<nucl16, std::char_traits<nucl16>> w{nucl16::A, nucl16::C, nucl16::G, nucl16::T, nucl16::U,
                                                          nucl16::N, nucl16::UNKNOWN};
    EXPECT_EQ(w, "ACGTUNN"_nucl16s);
}
