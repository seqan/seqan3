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
using structure_types  = ::testing::Types<dot_bracket3, wuss<>, dssp9>;

TYPED_TEST_CASE(structure, structure_types);

TYPED_TEST(structure, assign_char)
{
    using t = TypeParam;
    std::vector<char> input
    {
        '.', '(', ')',
        ':', ',', '-', '_', '~', ';',
        '<', '>', '[', ']', '{', '}',
        'H', 'B', 'E', 'G', 'I', 'T', 'S'
    };

    std::vector<TypeParam> cmp;
    if constexpr (std::is_same_v<TypeParam, dot_bracket3>)
    {
        cmp =
        {
            t::UNPAIRED, t::PAIR_OPEN, t::PAIR_CLOSE,
            t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN,
            t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN,
            t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN, t::UNKNOWN
        };
    }
    else if constexpr (std::is_same_v<TypeParam, wuss<>>)
    {
        cmp =
        {
            t::UNPAIRED, t::PAIR_OPEN1, t::PAIR_CLOSE1,
            t::UNPAIRED1, t::UNPAIRED2, t::UNPAIRED3, t::UNPAIRED4, t::UNPAIRED5, t::UNPAIRED6,
            t::PAIR_OPEN, t::PAIR_CLOSE, t::PAIR_OPEN2, t::PAIR_CLOSE2, t::PAIR_OPEN3, t::PAIR_CLOSE3,
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
            t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X, t::X,
            t::H, t::B, t::E, t::G, t::I, t::T, t::S
        };
    }

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
        EXPECT_EQ((assign_char(TypeParam{}, ch)), cm);
}

TYPED_TEST(structure, to_char)
{
    if constexpr (std::is_same_v<TypeParam, dot_bracket3>)
    {
        EXPECT_EQ(to_char(TypeParam::UNPAIRED), '.');
        EXPECT_EQ(to_char(TypeParam::PAIR_OPEN), '(');
        EXPECT_EQ(to_char(TypeParam::PAIR_CLOSE), ')');
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), '.');
    }
    else if constexpr (std::is_same_v<TypeParam, wuss<>>)
    {
        EXPECT_EQ(to_char(TypeParam::UNPAIRED), '.');
        EXPECT_EQ(to_char(TypeParam::UNPAIRED1), ':');
        EXPECT_EQ(to_char(TypeParam::UNPAIRED2), ',');
        EXPECT_EQ(to_char(TypeParam::UNPAIRED3), '-');
        EXPECT_EQ(to_char(TypeParam::UNPAIRED4), '_');
        EXPECT_EQ(to_char(TypeParam::UNPAIRED5), '~');
        EXPECT_EQ(to_char(TypeParam::UNPAIRED6), ';');
        EXPECT_EQ(to_char(TypeParam::PAIR_OPEN), '<');
        EXPECT_EQ(to_char(TypeParam::PAIR_CLOSE), '>');
        EXPECT_EQ(to_char(TypeParam::PAIR_OPEN1), '(');
        EXPECT_EQ(to_char(TypeParam::PAIR_CLOSE1), ')');
        EXPECT_EQ(to_char(TypeParam::PAIR_OPEN2), '[');
        EXPECT_EQ(to_char(TypeParam::PAIR_CLOSE2), ']');
        EXPECT_EQ(to_char(TypeParam::PAIR_OPEN3), '{');
        EXPECT_EQ(to_char(TypeParam::PAIR_CLOSE3), '}');
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

TEST(structure_stream_operator, dot_bracket3)
{
    std::stringstream ss;
    ss << dot_bracket3::PAIR_OPEN << dot_bracket3::PAIR_CLOSE << dot_bracket3::UNPAIRED;
    EXPECT_EQ(ss.str(), "().");
}

TEST(structure_stream_operator, wuss)
{
    std::stringstream ss;
    ss << wuss<>::PAIR_OPEN << wuss<>::PAIR_CLOSE << wuss<>::UNPAIRED;
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

TEST(dot_bracket3_literals, vector)
{
    std::vector<dot_bracket3> vec1;
    vec1.resize(5, dot_bracket3::PAIR_OPEN);
    EXPECT_EQ(vec1, "((((("_dot_bracket3);

    std::vector<dot_bracket3> vec2{dot_bracket3::UNPAIRED, dot_bracket3::PAIR_OPEN, dot_bracket3::PAIR_OPEN,
                                   dot_bracket3::PAIR_CLOSE, dot_bracket3::PAIR_CLOSE, dot_bracket3::UNKNOWN};
    EXPECT_EQ(vec2, ".(())."_dot_bracket3);
}

TEST(dot_bracket3_literals, basic_string)
{
    using string_t = std::basic_string<dot_bracket3, std::char_traits<dot_bracket3>>;
    string_t str1;
    str1.resize(5, dot_bracket3::PAIR_CLOSE);
    EXPECT_EQ(str1, ")))))"_dot_bracket3s);

    string_t str2{dot_bracket3::PAIR_OPEN, dot_bracket3::UNPAIRED, dot_bracket3::PAIR_CLOSE,
                  dot_bracket3::PAIR_OPEN, dot_bracket3::PAIR_CLOSE, dot_bracket3::UNKNOWN};
    EXPECT_EQ(str2, "(.)()."_dot_bracket3s);
}

TEST(wuss_literals, vector)
{
    std::vector<wuss<>> vec1;
    vec1.resize(5, wuss<>::PAIR_OPEN);
    EXPECT_EQ(vec1, "<<<<<"_wuss);

    std::vector<wuss<>> vec2{wuss<>::UNPAIRED, wuss<>::PAIR_OPEN, wuss<>::PAIR_OPEN,
                             wuss<>::PAIR_CLOSE, wuss<>::PAIR_CLOSE, wuss<>::UNPAIRED};
    EXPECT_EQ(vec2, ".<<>>."_wuss);
}

TEST(wuss_literals, basic_string)
{
    using string_t = std::basic_string<wuss<>, std::char_traits<wuss<>>>;
    string_t str1;
    str1.resize(5, wuss<>::PAIR_CLOSE);
    EXPECT_EQ(str1, ">>>>>"_wusss);

    string_t str2{wuss<>::PAIR_OPEN, wuss<>::UNPAIRED, wuss<>::PAIR_CLOSE,
                  wuss<>::PAIR_OPEN, wuss<>::PAIR_CLOSE, wuss<>::UNPAIRED};
    EXPECT_EQ(str2, "<.><>."_wusss);
}

TEST(dssp9_literals, vector)
{
    std::vector<dssp9> vec1;
    vec1.resize(5, dssp9::H);
    EXPECT_EQ(vec1, "HHHHH"_dssp9);

    std::vector<dssp9> vec2{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
    EXPECT_EQ(vec2, "EHHHTG"_dssp9);
}

TEST(dssp9_literals, basic_string)
{
    using string_t = std::basic_string<seqan3::dssp9, std::char_traits<seqan3::dssp9>>;
    string_t str1;
    str1.resize(5, dssp9::B);
    EXPECT_EQ(str1, "BBBBB"_dssp9s);

    string_t str2{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
    EXPECT_EQ(str2, "EHHHTG"_dssp9s);
}
