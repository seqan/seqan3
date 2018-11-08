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

template <typename T>
class nucleotide : public ::testing::Test
{};

// add all alphabets from the nucleotide sub module here
using nucleotide_types  = ::testing::Types<dna4, dna5, dna15, rna4, rna5, rna15>;
using nucleotide_types2 =       meta::list<dna4, dna5, dna15, rna4, rna5, rna15>; // needed for some tests

TYPED_TEST_CASE(nucleotide, nucleotide_types);

TYPED_TEST(nucleotide, assign_char_to_char)
{
    EXPECT_EQ(to_char(TypeParam{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(TypeParam{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(TypeParam{}.assign_char('G')), 'G');

    if constexpr (std::is_same_v<TypeParam, rna4> ||
                  std::is_same_v<TypeParam, rna5> ||
                  std::is_same_v<TypeParam, rna15>)
    {
        EXPECT_EQ(to_char(TypeParam{}.assign_char('U')), 'U');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('T')), 'U');
    }
    else if constexpr (std::is_same_v<TypeParam, dna4> ||
                       std::is_same_v<TypeParam, dna5> ||
                       std::is_same_v<TypeParam, dna15>)
    {
        EXPECT_EQ(to_char(TypeParam{}.assign_char('U')), 'T');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('T')), 'T');
    }

    if constexpr (alphabet_size_v<TypeParam> > 5)
    {
        EXPECT_EQ(to_char(TypeParam{}.assign_char('R')), 'R');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('Y')), 'Y');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('S')), 'S');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('W')), 'W');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('K')), 'K');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('M')), 'M');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('B')), 'B');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('D')), 'D');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('H')), 'H');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('V')), 'V');
    }

    if constexpr (std::is_same_v<TypeParam, dna4> || std::is_same_v<TypeParam, rna4>)
    {
        EXPECT_EQ(to_char(TypeParam{}.assign_char('N')), 'A');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('!')), 'A');
    }
    else
    {
        EXPECT_EQ(to_char(TypeParam{}.assign_char('N')), 'N');
        EXPECT_EQ(to_char(TypeParam{}.assign_char('!')), 'N');
    }
}

TYPED_TEST(nucleotide, complement)
{
    EXPECT_EQ(complement(TypeParam{}.assign_char('A')), TypeParam{}.assign_char('T'));
    EXPECT_EQ(complement(TypeParam{}.assign_char('C')), TypeParam{}.assign_char('G'));
    EXPECT_EQ(complement(TypeParam{}.assign_char('G')), TypeParam{}.assign_char('C'));
    EXPECT_EQ(complement(TypeParam{}.assign_char('T')), TypeParam{}.assign_char('A'));

    using vsize_t = std::decay_t<decltype(alphabet_size_v<TypeParam>)>;

    for (vsize_t i = 0u; i < alphabet_size_v<TypeParam>; ++i)
    {
        TypeParam c = assign_rank(TypeParam{}, i);

        EXPECT_EQ(complement(complement(c)), c);
    }
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
    EXPECT_EQ(other_type{TypeParam{}.assign_char('C')}, other_type{}.assign_char('C'));
    // assign
    other_type l{};
    l = TypeParam{}.assign_char('C');
    EXPECT_EQ(l, other_type{}.assign_char('C'));
}

// conversion to any other nucleotide type
TYPED_TEST(nucleotide, explicit_conversion)
{
    meta::for_each(nucleotide_types2{}, [&] (auto && nucl) constexpr
    {
        using out_type = std::decay_t<decltype(nucl)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('A')), out_type{}.assign_char('A'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('C')), out_type{}.assign_char('C'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('G')), out_type{}.assign_char('G'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('T'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('U')), out_type{}.assign_char('U'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('U'));
    });
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(dna4_literals, vector)
{
    dna4_vector v;
    v.resize(5, 'A'_dna4);
    EXPECT_EQ(v, "AAAAA"_dna4);

    std::vector<dna4> w{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'U'_dna4, 'N'_dna4};
    EXPECT_EQ(w, "ACGTTA"_dna4);
}

TEST(dna5_literals, vector)
{
    dna5_vector v;
    v.resize(5, 'A'_dna5);
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<dna5> w{'A'_dna5, 'C'_dna5, 'G'_dna5, 'T'_dna5, 'U'_dna5, 'N'_dna5, 'N'_dna5};
    EXPECT_EQ(w, "ACGTTNN"_dna5);
}

TEST(dna15_literals, vector)
{
    dna15_vector v;
    v.resize(5, 'A'_dna15);
    EXPECT_EQ(v, "AAAAA"_dna15);

    std::vector<dna15> w{'A'_dna15, 'C'_dna15, 'G'_dna15, 'T'_dna15, 'U'_dna15, 'N'_dna15, 'N'_dna15};
    EXPECT_EQ(w, "ACGTTNN"_dna15);
}

TEST(rna4_literals, vector)
{
    rna4_vector v;
    v.resize(5, 'A'_rna4);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<rna4> w{'A'_rna4, 'C'_rna4, 'G'_rna4, 'T'_rna4, 'U'_rna4, 'N'_rna4};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

TEST(rna5_literals, vector)
{
    rna5_vector v;
    v.resize(5, 'A'_rna5);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<rna5> w{'A'_rna5, 'C'_rna5, 'G'_rna5, 'T'_rna5, 'U'_rna5, 'N'_rna5, 'N'_rna5};
    EXPECT_EQ(w, "ACGUUNN"_rna5);
}

TEST(rna15_literals, vector)
{
    rna15_vector v;
    v.resize(5, 'A'_rna15);
    EXPECT_EQ(v, "AAAAA"_rna15);

    std::vector<rna15> w{'A'_rna15, 'C'_rna15, 'G'_rna15, 'T'_rna15, 'U'_rna15, 'N'_rna15, 'N'_rna15};
    EXPECT_EQ(w, "ACGTTNN"_rna15);
}
