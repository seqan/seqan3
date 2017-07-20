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

/*!\file
 * \ingroup alphabet
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Tests for the structure alphabets.
 */

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/structure/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class structure : public ::testing::Test
{};

// add all alphabets from the structure sub module here
using structure_types  = ::testing::Types<db3, wuss<>>;

TYPED_TEST_CASE(structure, structure_types);

TYPED_TEST(structure, assign_char)
{
    using t = TypeParam;
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '_', '~', '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<TypeParam> cmp;
    if constexpr (std::is_same_v<TypeParam, db3>)
    {
        cmp =
        {
            t::NP, t::BL, t::BR,
            t::NA, t::NA, t::NA, t::NA, t::NA, t::NA, t::NA, t::NA, t::NA, t::NA,
            t::NA, t::NA, t::NA, t::NA, t::NA, t::NA, t::NA
        };
    }
    else if constexpr (std::is_same_v<TypeParam, wuss<>>)
    {
        cmp =
        {
            t::NP, t::BL1, t::BR1,
            t::NP1, t::NP2, t::NP3, t::NP4, t::BL, t::BR, t::BL2, t::BR2, t::BL3, t::BR3,
            "H"_wuss.front(),
            "B"_wuss.front(),
            "E"_wuss.front(),
            "G"_wuss.front(),
            "I"_wuss.front(),
            "T"_wuss.front(),
            "S"_wuss.front()
        };
    }
    else if constexpr (std::is_same_v<TypeParam, dssp9>)
    {
        cmp =
        {
            t::X, t::X, t::X,
            t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X,
            t::H, t::B, t::E, t::G, t::I, t::T, t::S
        };
    }

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
    EXPECT_EQ((assign_char(TypeParam{}, ch)), cm);
}

TYPED_TEST(structure, to_char)
{
    if constexpr (std::is_same_v<TypeParam, db3>)
    {
        EXPECT_EQ(to_char(TypeParam::NP), '.');
        EXPECT_EQ(to_char(TypeParam::BL), '(');
        EXPECT_EQ(to_char(TypeParam::BR), ')');
        EXPECT_EQ(to_char(TypeParam::NA), '.');
    }
    else if constexpr (std::is_same_v<TypeParam, wuss<>>)
    {
        EXPECT_EQ(to_char(TypeParam::NP), '.');
        EXPECT_EQ(to_char(TypeParam::NP1), ':');
        EXPECT_EQ(to_char(TypeParam::NP2), ',');
        EXPECT_EQ(to_char(TypeParam::NP3), '_');
        EXPECT_EQ(to_char(TypeParam::NP4), '~');
        EXPECT_EQ(to_char(TypeParam::BL), '<');
        EXPECT_EQ(to_char(TypeParam::BR), '>');
        EXPECT_EQ(to_char(TypeParam::BL1), '(');
        EXPECT_EQ(to_char(TypeParam::BR1), ')');
        EXPECT_EQ(to_char(TypeParam::BL2), '[');
        EXPECT_EQ(to_char(TypeParam::BR2), ']');
        EXPECT_EQ(to_char(TypeParam::BL3), '{');
        EXPECT_EQ(to_char(TypeParam::BR3), '}');
    }
    else if constexpr (std::is_same_v<TypeParam, dssp9>)
    {
        EXPECT_EQ(to_char(TypeParam::H), 'H');
        EXPECT_EQ(to_char(TypeParam::B), 'B');
        EXPECT_EQ(to_char(TypeParam::E), 'E');
        EXPECT_EQ(to_char(TypeParam::G), 'G');
        EXPECT_EQ(to_char(TypeParam::I), 'I');
        EXPECT_EQ(to_char(TypeParam::T), 'T');
        EXPECT_EQ(to_char(TypeParam::S), 'S');
        EXPECT_EQ(to_char(TypeParam::C), 'C');
        EXPECT_EQ(to_char(TypeParam::X), 'X');
    }
}

TYPED_TEST(structure, concept)
{
    EXPECT_TRUE(structure_concept<TypeParam>);
}

TEST(structure_stream_operator, db3)
{
    std::stringstream ss;
    ss << db3::BL << db3::BR << db3::NP;
    EXPECT_EQ(ss.str(), "().");
}

TEST(structure_stream_operator, wuss)
{
    std::stringstream ss;
    ss << wuss<>::BL << wuss<>::BR << wuss<>::NP;
    EXPECT_EQ(ss.str(), "<>.");
}

TEST(structure_stream_operator, dssp9)
{
    std::stringstream ss;
    ss << dssp9::E << dssp9::H << dssp9::C;
    EXPECT_EQ(ss.str(), "EHC");
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(db3_literals, vector)
{
    db3_vector v;
    v.resize(5, db3::BL);
    EXPECT_EQ(v, "((((("_db3);

    std::vector<db3> w{db3::NP, db3::BL, db3::BL, db3::BR, db3::BR, db3::NA};
    EXPECT_EQ(w, ".(())."_db3);
}

TEST(db3_literals, basic_string)
{
    db3_string v;
    v.resize(5, db3::BR);
    EXPECT_EQ(v, ")))))"_db3s);

    std::basic_string<db3, std::char_traits<db3>> w{db3::BL, db3::NP, db3::BR, db3::BL, db3::BR, db3::NA};
    EXPECT_EQ(w, "(.)()."_db3s);
}

TEST(wuss_literals, vector)
{
    wuss_vector v;
    v.resize(5, wuss<>::BL);
    EXPECT_EQ(v, "<<<<<"_wuss);

    std::vector<wuss<>> w{wuss<>::NP, wuss<>::BL, wuss<>::BL, wuss<>::BR, wuss<>::BR, wuss<>::NP};
    EXPECT_EQ(w, ".<<>>."_wuss);
}

TEST(wuss_literals, basic_string)
{
    wuss_string v;
    v.resize(5, wuss<>::BR);
    EXPECT_EQ(v, ">>>>>"_wusss);

    std::basic_string<wuss<>, std::char_traits<wuss<>>> w{wuss<>::BL, wuss<>::NP, wuss<>::BR, wuss<>::BL, wuss<>::BR,
                                                          wuss<>::NP};
    EXPECT_EQ(w, "<.><>."_wusss);
}

TEST(dssp9_literals, vector)
{
    dssp9_vector v;
    v.resize(5, dssp9::H);
    EXPECT_EQ(v, "HHHHH"_dssp9);

    std::vector<dssp9> w{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
    EXPECT_EQ(w, "EHHHTG"_dssp9);
}

TEST(dssp9_literals, basic_string)
{
    dssp9_string v;
    v.resize(5, dssp9::B);
    EXPECT_EQ(v, "BBBBB"_dssp9s);

    std::basic_string<dssp9, std::char_traits<dssp9>> w{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
    EXPECT_EQ(w, "EHHHTG"_dssp9s);
}
