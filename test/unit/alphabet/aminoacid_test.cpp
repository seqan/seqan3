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
using aminoacid_types = ::testing::Types<aa27, aa20>;
using aminoacid_types2 =      meta::list<aa27, aa20>; // needed for some tests

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
    std::vector<TypeParam> cmp;
    if constexpr (std::is_same_v<TypeParam, aa27>)
    {
        cmp =
        {
            t::A, t::B, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::J, t::K, t::L, t::M,
            t::A, t::B, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::J, t::K, t::L, t::M,
            t::N, t::O, t::P, t::Q, t::R, t::S, t::T, t::U, t::V, t::W, t::X, t::Y, t::Z,
            t::N, t::O, t::P, t::Q, t::R, t::S, t::T, t::U, t::V, t::W, t::X, t::Y, t::Z,
            t::TERMINATOR, t::X
        };
    }
    else if constexpr (std::is_same_v<TypeParam, aa20>)
    {
        cmp =
        {
            t::A, t::D, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::L, t::K, t::L, t::M,
            t::A, t::D, t::C, t::D, t::E, t::F, t::G, t::H, t::I, t::L, t::K, t::L, t::M,
            t::N, t::L, t::P, t::Q, t::R, t::S, t::T, t::C, t::V, t::W, t::S, t::Y, t::E,
            t::N, t::L, t::P, t::Q, t::R, t::S, t::T, t::C, t::V, t::W, t::S, t::Y, t::E,
            t::W, t::S
        };
    }

    for (auto [ ch, cm ] : ranges::view::zip(input, cmp))
        EXPECT_EQ((assign_char(TypeParam{}, ch)), cm);
}

TYPED_TEST(aminoacid, to_char)
{
    EXPECT_EQ(to_char(TypeParam::A), 'A');
    EXPECT_EQ(to_char(TypeParam::C), 'C');
    EXPECT_EQ(to_char(TypeParam::D), 'D');
    EXPECT_EQ(to_char(TypeParam::E), 'E');
    EXPECT_EQ(to_char(TypeParam::F), 'F');
    EXPECT_EQ(to_char(TypeParam::G), 'G');
    EXPECT_EQ(to_char(TypeParam::H), 'H');
    EXPECT_EQ(to_char(TypeParam::I), 'I');
    EXPECT_EQ(to_char(TypeParam::K), 'K');
    EXPECT_EQ(to_char(TypeParam::L), 'L');
    EXPECT_EQ(to_char(TypeParam::M), 'M');
    EXPECT_EQ(to_char(TypeParam::N), 'N');
    EXPECT_EQ(to_char(TypeParam::P), 'P');
    EXPECT_EQ(to_char(TypeParam::Q), 'Q');
    EXPECT_EQ(to_char(TypeParam::R), 'R');
    EXPECT_EQ(to_char(TypeParam::S), 'S');
    EXPECT_EQ(to_char(TypeParam::T), 'T');
    EXPECT_EQ(to_char(TypeParam::V), 'V');
    EXPECT_EQ(to_char(TypeParam::W), 'W');
    EXPECT_EQ(to_char(TypeParam::Y), 'Y');

    if constexpr (std::is_same_v<TypeParam, aa27>)
    {
        EXPECT_EQ(to_char(TypeParam::B), 'B');
        EXPECT_EQ(to_char(TypeParam::J), 'J');
        EXPECT_EQ(to_char(TypeParam::O), 'O');
        EXPECT_EQ(to_char(TypeParam::U), 'U');
        EXPECT_EQ(to_char(TypeParam::X), 'X');
        EXPECT_EQ(to_char(TypeParam::Z), 'Z');
        EXPECT_EQ(to_char(TypeParam::TERMINATOR), '*');
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'X');
    }
    else if constexpr (std::is_same_v<TypeParam, aa20>)
    {
        EXPECT_EQ(to_char(TypeParam::B), 'D');
        EXPECT_EQ(to_char(TypeParam::J), 'L');
        EXPECT_EQ(to_char(TypeParam::O), 'L');
        EXPECT_EQ(to_char(TypeParam::U), 'C');
        EXPECT_EQ(to_char(TypeParam::X), 'S');
        EXPECT_EQ(to_char(TypeParam::Z), 'E');
        EXPECT_EQ(to_char(TypeParam::TERMINATOR), 'W');
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), 'S');
    }
}

TYPED_TEST(aminoacid, stream_operator)
{
    std::stringstream ss;
    ss << TypeParam::A << TypeParam::C << TypeParam::G << TypeParam::B << TypeParam::J
       << TypeParam::O << TypeParam::U << TypeParam::X << TypeParam::Z;

    if constexpr (std::is_same_v<TypeParam, aa27>)
    {
        EXPECT_EQ(ss.str(), "ACGBJOUXZ");
    }
    else if constexpr (std::is_same_v<TypeParam, aa20>)
    {
        EXPECT_EQ(ss.str(), "ACGDLLCSE");
    }
}

TYPED_TEST(aminoacid, concept_check)
{
    EXPECT_TRUE(aminoacid_concept<TypeParam>);
    EXPECT_TRUE(aminoacid_concept<TypeParam &>);
}

// ------------------------------------------------------------------
// conversions
// ------------------------------------------------------------------

// conversion to any other amino acid type
TYPED_TEST(aminoacid, explicit_conversion)
{
    meta::for_each(aminoacid_types2{}, [&] (auto && aa) constexpr
    {
        using out_type = std::decay_t<decltype(aa)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam::A), out_type::A);
        EXPECT_EQ(static_cast<out_type>(TypeParam::C), out_type::C);
        EXPECT_EQ(static_cast<out_type>(TypeParam::D), out_type::D);
        EXPECT_EQ(static_cast<out_type>(TypeParam::E), out_type::E);
        EXPECT_EQ(static_cast<out_type>(TypeParam::F), out_type::F);
        EXPECT_EQ(static_cast<out_type>(TypeParam::G), out_type::G);
        EXPECT_EQ(static_cast<out_type>(TypeParam::H), out_type::H);
        EXPECT_EQ(static_cast<out_type>(TypeParam::I), out_type::I);
        EXPECT_EQ(static_cast<out_type>(TypeParam::K), out_type::K);
        EXPECT_EQ(static_cast<out_type>(TypeParam::L), out_type::L);
        EXPECT_EQ(static_cast<out_type>(TypeParam::M), out_type::M);
        EXPECT_EQ(static_cast<out_type>(TypeParam::N), out_type::N);
        EXPECT_EQ(static_cast<out_type>(TypeParam::P), out_type::P);
        EXPECT_EQ(static_cast<out_type>(TypeParam::Q), out_type::Q);
        EXPECT_EQ(static_cast<out_type>(TypeParam::R), out_type::R);
        EXPECT_EQ(static_cast<out_type>(TypeParam::S), out_type::S);
        EXPECT_EQ(static_cast<out_type>(TypeParam::T), out_type::T);
        EXPECT_EQ(static_cast<out_type>(TypeParam::V), out_type::V);
        EXPECT_EQ(static_cast<out_type>(TypeParam::W), out_type::W);
        EXPECT_EQ(static_cast<out_type>(TypeParam::Y), out_type::Y);
        if (std::is_same_v<TypeParam, aa27>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam::B), out_type::B);
            EXPECT_EQ(static_cast<out_type>(TypeParam::J), out_type::J);
            EXPECT_EQ(static_cast<out_type>(TypeParam::O), out_type::O);
            EXPECT_EQ(static_cast<out_type>(TypeParam::U), out_type::U);
            EXPECT_EQ(static_cast<out_type>(TypeParam::X), out_type::X);
            EXPECT_EQ(static_cast<out_type>(TypeParam::Z), out_type::Z);
            EXPECT_EQ(static_cast<out_type>(TypeParam::TERMINATOR), out_type::TERMINATOR);
            EXPECT_EQ(static_cast<out_type>(TypeParam::UNKNOWN), out_type::UNKNOWN);
        }
        else if (std::is_same_v<TypeParam, aa20>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam::B), out_type::D);
            EXPECT_EQ(static_cast<out_type>(TypeParam::J), out_type::L);
            EXPECT_EQ(static_cast<out_type>(TypeParam::O), out_type::L);
            EXPECT_EQ(static_cast<out_type>(TypeParam::U), out_type::C);
            EXPECT_EQ(static_cast<out_type>(TypeParam::X), out_type::S);
            EXPECT_EQ(static_cast<out_type>(TypeParam::Z), out_type::E);
            EXPECT_EQ(static_cast<out_type>(TypeParam::TERMINATOR), out_type::W);
            EXPECT_EQ(static_cast<out_type>(TypeParam::UNKNOWN), out_type::S);
        }
    });
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(literals, vector)
{
    aa27_vector v27;
    aa20_vector v20;
    v27.resize(5, aa27::A);
    v20.resize(5, aa20::B);
    EXPECT_EQ(v27, "AAAAA"_aa27);
    EXPECT_EQ(v20, "DDDDD"_aa20);

    std::vector<aa27> w27{aa27::A, aa27::Y, aa27::P, aa27::T, aa27::U, aa27::N, aa27::X, aa27::UNKNOWN,
                        aa27::TERMINATOR};
    std::vector<aa20> w20{aa20::A, aa20::B, aa20::J, aa20::O, aa20::U, aa20::X, aa20::Z, aa20::UNKNOWN,
                        aa20::TERMINATOR, aa20::TERMINATOR};
    EXPECT_EQ(w27, "AYPTUNXX*"_aa27);
    EXPECT_EQ(w20, "ADLLCSESW*"_aa20);
}

// ------------------------------------------------------------------
// comparators
// ------------------------------------------------------------------

TYPED_TEST(aminoacid, comparators)
{
    EXPECT_TRUE(TypeParam::A == TypeParam::A);
    EXPECT_TRUE(TypeParam::A != TypeParam::B);
    EXPECT_TRUE(TypeParam::A < TypeParam::B);
    EXPECT_TRUE(TypeParam::A <= TypeParam::B);
    EXPECT_TRUE(TypeParam::B > TypeParam::A);
    EXPECT_TRUE(TypeParam::B >= TypeParam::A);
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
    auto range_triplet = ranges::view::concat(ranges::view::single(n1), ranges::view::single(n2),
                                              ranges::view::single(n3));
    aa27 t2{translate_triplet(range_triplet)};

    EXPECT_EQ(t2, c);

    // Tuple interface
    std::tuple tuple_triplet{n1, n2, n3};
    aa27 t3{translate_triplet(tuple_triplet)};

    EXPECT_EQ(t3, c);
}
