// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class nucleotide : public ::testing::Test
{};

// add all alphabets from the nucleotide sub module here
using nucleotide_types  = ::testing::Types<dna4, dna5, dna15, rna4, rna5, rna15>;
using nucleotide_types2 =       meta::list<dna4, dna5, dna15, rna4, rna5, rna15>; // needed for some tests

TYPED_TEST_CASE(nucleotide, nucleotide_types);

TYPED_TEST(nucleotide, assign_char)
{
    using t = TypeParam;
    std::vector<char> input
    {
        'A', 'C', 'G', 'T', 'U', 'N',
        'a', 'c', 'g', 't', 'u', 'n',
        'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V',
        'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
        '!'
    };

    std::vector<TypeParam> cmp;
    if constexpr (std::is_same_v<TypeParam, dna4> || std::is_same_v<TypeParam, rna4>)
    {
        cmp =
        {
            t::A, t::C, t::G, t::T, t::U, t::A,
            t::A, t::C, t::G, t::T, t::U, t::A,
            t::A, t::C, t::C, t::A, t::G, t::A, t::C, t::A, t::A, t::A,
            t::A, t::C, t::C, t::A, t::G, t::A, t::C, t::A, t::A, t::A,
            t::A
        };
    } else if constexpr (std::is_same_v<TypeParam, dna5> || std::is_same_v<TypeParam, rna5>)
    {
        cmp =
        {
            t::A, t::C, t::G, t::T, t::U, t::N,
            t::A, t::C, t::G, t::T, t::U, t::N,
            t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N,
            t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N, t::N,
            t::N
        };
    } else if constexpr (std::is_same_v<TypeParam, dna15> || std::is_same_v<TypeParam, rna15>)
    {
        cmp =
        {
            t::A, t::C, t::G, t::T, t::U, t::N,
            t::A, t::C, t::G, t::T, t::U, t::N,
            t::R, t::Y, t::S, t::W, t::K, t::M, t::B, t::D, t::H, t::V,
            t::R, t::Y, t::S, t::W, t::K, t::M, t::B, t::D, t::H, t::V,
            t::N
        };
    }

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
        EXPECT_EQ((assign_char(TypeParam{}, ch)), cm);
}

TYPED_TEST(nucleotide, to_char)
{
    EXPECT_EQ(to_char(TypeParam::A), 'A');
    EXPECT_EQ(to_char(TypeParam::C), 'C');
    EXPECT_EQ(to_char(TypeParam::G), 'G');

    if constexpr (std::is_same_v<TypeParam, rna4> ||
                  std::is_same_v<TypeParam, rna5> ||
                  std::is_same_v<TypeParam, rna15>)
    {
        EXPECT_EQ(to_char(TypeParam::U), 'U');
        EXPECT_EQ(to_char(TypeParam::T), 'U');
    }
    else if constexpr (std::is_same_v<TypeParam, dna4> ||
                       std::is_same_v<TypeParam, dna5> ||
                       std::is_same_v<TypeParam, dna15>)
    {
        EXPECT_EQ(to_char(TypeParam::U), 'T');
        EXPECT_EQ(to_char(TypeParam::T), 'T');
    }

    if constexpr (alphabet_size_v<TypeParam> > 5)
    {
        EXPECT_EQ(to_char(TypeParam::R), 'R');
        EXPECT_EQ(to_char(TypeParam::Y), 'Y');
        EXPECT_EQ(to_char(TypeParam::S), 'S');
        EXPECT_EQ(to_char(TypeParam::W), 'W');
        EXPECT_EQ(to_char(TypeParam::K), 'K');
        EXPECT_EQ(to_char(TypeParam::M), 'M');
        EXPECT_EQ(to_char(TypeParam::B), 'B');
        EXPECT_EQ(to_char(TypeParam::D), 'D');
        EXPECT_EQ(to_char(TypeParam::H), 'H');
        EXPECT_EQ(to_char(TypeParam::V), 'V');
    }

    if constexpr (std::is_same_v<TypeParam, dna4> || std::is_same_v<TypeParam, rna4>)
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'A');
    else
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'N');
}

TYPED_TEST(nucleotide, complement)
{
    EXPECT_EQ(complement(TypeParam::A), TypeParam::T);
    EXPECT_EQ(complement(TypeParam::C), TypeParam::G);
    EXPECT_EQ(complement(TypeParam::G), TypeParam::C);
    EXPECT_EQ(complement(TypeParam::T), TypeParam::A);

    using vsize_t = std::decay_t<decltype(alphabet_size_v<TypeParam>)>;

    for (vsize_t i = 0u; i < alphabet_size_v<TypeParam>; ++i)
    {
        TypeParam c = assign_rank(TypeParam{}, i);

        EXPECT_EQ(complement(complement(c)), c);
    }
}

TYPED_TEST(nucleotide, stream_operator)
{
    std::stringstream ss;
    ss << TypeParam::A << TypeParam::C << TypeParam::G;
    EXPECT_EQ(ss.str(), "ACG");
}

TYPED_TEST(nucleotide, concept_check)
{
    EXPECT_TRUE(nucleotide_concept<TypeParam>);
    EXPECT_TRUE(nucleotide_concept<TypeParam &>);
}

// ------------------------------------------------------------------
// conversion
// ------------------------------------------------------------------

// conversion to rna/dna of same size
TYPED_TEST(nucleotide, implicit_conversion)
{
    using other_type = std::conditional_t<std::is_same_v<TypeParam, rna4>,  dna4,
                       std::conditional_t<std::is_same_v<TypeParam, dna4>,  rna4,
                       std::conditional_t<std::is_same_v<TypeParam, rna5>,  dna5,
                       std::conditional_t<std::is_same_v<TypeParam, dna5>,  rna5,
                       std::conditional_t<std::is_same_v<TypeParam, dna15>, rna15,
                       /* must be rna15 */                                  dna15>>>>>;

    // construct
    EXPECT_EQ(other_type{TypeParam::C}, other_type::C);
    // assign
    other_type l{};
    l = TypeParam::C;
    EXPECT_EQ(l, other_type::C);
}

// conversion to any other nucleotide type
TYPED_TEST(nucleotide, explicit_conversion)
{
    meta::for_each(nucleotide_types2{}, [&] (auto && nucl) constexpr
    {
        using out_type = std::decay_t<decltype(nucl)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam::A), out_type::A);
        EXPECT_EQ(static_cast<out_type>(TypeParam::C), out_type::C);
        EXPECT_EQ(static_cast<out_type>(TypeParam::G), out_type::G);
        EXPECT_EQ(static_cast<out_type>(TypeParam::T), out_type::T);
        EXPECT_EQ(static_cast<out_type>(TypeParam::U), out_type::U);
        EXPECT_EQ(static_cast<out_type>(TypeParam::T), out_type::U);
    });
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

TEST(dna5_literals, vector)
{
    dna5_vector v;
    v.resize(5, dna5::A);
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<dna5> w{dna5::A, dna5::C, dna5::G, dna5::T, dna5::U, dna5::N, dna5::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_dna5);
}

TEST(dna15_literals, vector)
{
    dna15_vector v;
    v.resize(5, dna15::A);
    EXPECT_EQ(v, "AAAAA"_dna15);

    std::vector<dna15> w{dna15::A, dna15::C, dna15::G, dna15::T, dna15::U, dna15::N, dna15::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_dna15);
}

TEST(rna4_literals, vector)
{
    rna4_vector v;
    v.resize(5, rna4::A);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<rna4> w{rna4::A, rna4::C, rna4::G, rna4::T, rna4::U, rna4::UNKNOWN};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

TEST(rna5_literals, vector)
{
    rna5_vector v;
    v.resize(5, rna5::A);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<rna5> w{rna5::A, rna5::C, rna5::G, rna5::T, rna5::U, rna5::N, rna5::UNKNOWN};
    EXPECT_EQ(w, "ACGUUNN"_rna5);
}

TEST(rna15_literals, vector)
{
    rna15_vector v;
    v.resize(5, rna15::A);
    EXPECT_EQ(v, "AAAAA"_rna15);

    std::vector<rna15> w{rna15::A, rna15::C, rna15::G, rna15::T, rna15::U, rna15::N, rna15::UNKNOWN};
    EXPECT_EQ(w, "ACGTTNN"_rna15);
}
