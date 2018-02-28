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

#include <range/v3/view/concat.hpp>
#include <range/v3/view/single.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class aminoacid : public ::testing::Test
{};

// add all alphabets from the aminoacid sub module here
using aminoacid_types  = ::testing::Types<aa27>;

TYPED_TEST_CASE(aminoacid, aminoacid_types);

TYPED_TEST(aminoacid, assign_char)
{
    using t = TypeParam;
    std::vector<char> input
    {
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        '*', '!'
    };

    std::vector<TypeParam> cmp
    {
        t::A, t::B, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::J, t::K, t::L, t::M,
        t::A, t::B, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::J, t::K, t::L, t::M,
        t::N, t::O, t::P, t::Q, t::R, t::S, t::T, t::U, t::V, t::W, t::X, t::Y, t::Z,
        t::N, t::O, t::P, t::Q, t::R, t::S, t::T, t::U, t::V, t::W, t::X, t::Y, t::Z,
        t::TERMINATOR, t::X
    };

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
        EXPECT_EQ((assign_char(TypeParam{}, ch)), cm);
}

TYPED_TEST(aminoacid, to_char)
{
    EXPECT_EQ(to_char(TypeParam::A), 'A');
    EXPECT_EQ(to_char(TypeParam::B), 'B');
    EXPECT_EQ(to_char(TypeParam::C), 'C');
    EXPECT_EQ(to_char(TypeParam::D), 'D');
    EXPECT_EQ(to_char(TypeParam::E), 'E');
    EXPECT_EQ(to_char(TypeParam::F), 'F');
    EXPECT_EQ(to_char(TypeParam::G), 'G');
    EXPECT_EQ(to_char(TypeParam::H), 'H');
    EXPECT_EQ(to_char(TypeParam::I), 'I');
    EXPECT_EQ(to_char(TypeParam::J), 'J');
    EXPECT_EQ(to_char(TypeParam::K), 'K');
    EXPECT_EQ(to_char(TypeParam::L), 'L');
    EXPECT_EQ(to_char(TypeParam::M), 'M');
    EXPECT_EQ(to_char(TypeParam::N), 'N');
    EXPECT_EQ(to_char(TypeParam::O), 'O');
    EXPECT_EQ(to_char(TypeParam::P), 'P');
    EXPECT_EQ(to_char(TypeParam::Q), 'Q');
    EXPECT_EQ(to_char(TypeParam::R), 'R');
    EXPECT_EQ(to_char(TypeParam::S), 'S');
    EXPECT_EQ(to_char(TypeParam::T), 'T');
    EXPECT_EQ(to_char(TypeParam::U), 'U');
    EXPECT_EQ(to_char(TypeParam::V), 'V');
    EXPECT_EQ(to_char(TypeParam::W), 'W');
    EXPECT_EQ(to_char(TypeParam::X), 'X');
    EXPECT_EQ(to_char(TypeParam::Y), 'Y');
    EXPECT_EQ(to_char(TypeParam::Z), 'Z');
    EXPECT_EQ(to_char(TypeParam::TERMINATOR), '*');
    EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'X');
}

TYPED_TEST(aminoacid, stream_operator)
{
    std::stringstream ss;
    ss << TypeParam::A << TypeParam::C << TypeParam::G;
    EXPECT_EQ(ss.str(), "ACG");
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(aa27_literals, vector)
{
    aa27_vector v;
    v.resize(5, aa27::A);
    EXPECT_EQ(v, "AAAAA"_aa27);

    std::vector<aa27> w{aa27::A, aa27::Y, aa27::P, aa27::T, aa27::U, aa27::N, aa27::X, aa27::UNKNOWN,
                        aa27::TERMINATOR};
    EXPECT_EQ(w, "AYPTUNXX*"_aa27);
}

TEST(aa27_literals, string)
{
    aa27_string v;
    v.resize(5, aa27::A);
    EXPECT_EQ(v, "AAAAA"_aa27s);

    std::basic_string<aa27, std::char_traits<aa27>> w{aa27::A, aa27::Y, aa27::P, aa27::T, aa27::U,
                                                      aa27::N, aa27::X, aa27::UNKNOWN, aa27::TERMINATOR};
    EXPECT_EQ(w, "AYPTUNXX*"_aa27s);
}

// ------------------------------------------------------------------
// translation
// ------------------------------------------------------------------

TEST(translation, translate_triplets)
{
    dna15 n1{dna15::C};
    dna15 n2{dna15::T};
    dna15 n3{dna15::A};
    aa27 c{aa27::L};

    // Nucleotide interface
    aa27 t1{translate_triplet<genetic_code::CANONICAL, dna15>(n1, n2, n3)};

    EXPECT_EQ(t1, c);

    // Range interface
    auto range_triplet = ranges::view::concat(ranges::view::single(n1), ranges::view::single(n2), ranges::view::single(n3));
    aa27 t2{translate_triplet(range_triplet)};

    EXPECT_EQ(t2, c);

    // Tuple interface
    std::tuple tuple_triplet{n1, n2, n3};
    aa27 t3{translate_triplet(tuple_triplet)};

    EXPECT_EQ(t3, c);
}
