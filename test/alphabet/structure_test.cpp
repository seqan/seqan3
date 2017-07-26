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

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
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

// ------------------------------------------------------------------
// composition nucleotide x structure
// ------------------------------------------------------------------

// default/zero construction
TEST(structure_composition, ctr)
{
    [[maybe_unused]] structure_composition<rna4, dot_bracket3> t1;
}

// aggregate initialization
TEST(structure_composition, aggr)
{
    structure_composition<rna4, dot_bracket3> t1;
    structure_composition<rna4, dot_bracket3> t2{rna4::C, dot_bracket3::PAIR_CLOSE};
    EXPECT_NE(t1, t2);
}

// zero initialization
TEST(structure_composition, zro)
{
    structure_composition<rna4, dot_bracket3> t1{rna4::A, dot_bracket3::UNKNOWN};
    structure_composition<rna4, dot_bracket3> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(structure_composition, cp_ctr)
{
    structure_composition<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_OPEN};
    structure_composition<rna4, dot_bracket3> t2{t1};
    structure_composition<rna4, dot_bracket3> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(structure_composition, mv_ctr)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::UNPAIRED};
    structure_composition<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::UNPAIRED};
    structure_composition<rna4, dot_bracket3> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    structure_composition<rna4, dot_bracket3> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(structure_composition, cp_assgn)
{
    structure_composition<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::UNPAIRED};
    structure_composition<rna4, dot_bracket3> t2;
    structure_composition<rna4, dot_bracket3> t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(structure_composition, mv_assgn)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::UNPAIRED};
    structure_composition<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::UNPAIRED};
    structure_composition<rna4, dot_bracket3> t2;
    structure_composition<rna4, dot_bracket3> t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(structure_composition, swap)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};
    structure_composition<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_OPEN};
    structure_composition<rna4, dot_bracket3> t2{};
    structure_composition<rna4, dot_bracket3> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// get<1>
TEST(structure_composition, get_i)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), rna4 &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), dot_bracket3 &>);

    EXPECT_EQ(seqan3::get<0>(t0), rna4::C);
    EXPECT_EQ(seqan3::get<1>(t0), dot_bracket3{dot_bracket3::PAIR_OPEN});
}

// std::get<1>
TEST(structure_composition, stdget_i)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::UNPAIRED};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), rna4 &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), dot_bracket3 &>);

    EXPECT_EQ(std::get<0>(t0), rna4::C);
    EXPECT_EQ(std::get<1>(t0), dot_bracket3{dot_bracket3::UNPAIRED});
}

// structured bindings
TEST(structure_composition, struct_binding)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};
    auto [ i, l ] = t0;

    static_assert(std::is_same_v<decltype(i), rna4>);
    static_assert(std::is_same_v<decltype(l), dot_bracket3>);

    EXPECT_EQ(i, rna4::C);
    EXPECT_EQ(l, dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// get<type>
TEST(structure_composition, get_type)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};

    EXPECT_EQ(seqan3::get<rna4>(t0), rna4::C);
    EXPECT_EQ(seqan3::get<dot_bracket3>(t0), dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// std::get<type>
TEST(structure_composition, stdget_type)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};

    EXPECT_EQ(std::get<rna4>(t0), rna4::C);
    EXPECT_EQ(std::get<dot_bracket3>(t0), dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// std::tuple_element
TEST(structure_composition, tuple_element)
{
    using pt = structure_composition<rna4, dot_bracket3>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, rna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, dot_bracket3>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// type deduction
TEST(structure_composition, type_deduce)
{
    structure_composition t0{rna4::C, dot_bracket3{dot_bracket3::PAIR_CLOSE}};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, rna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, dot_bracket3>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// explicit cast to element
TEST(structure_composition, cast_to_element)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};

    auto d = static_cast<rna4>(t0);
    auto q = static_cast<dot_bracket3>(t0);
    static_assert(std::is_same_v<decltype(d), rna4>);
    static_assert(std::is_same_v<decltype(q), dot_bracket3>);

    EXPECT_EQ(d, rna4::C);
    EXPECT_EQ(q, dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// comparison operators
TEST(structure_composition, cmp)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};
    structure_composition<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_CLOSE};
    structure_composition<rna4, dot_bracket3> t2{rna4::G, dot_bracket3::PAIR_CLOSE};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

// alphabet_concept: rank_type
TEST(structure_composition, rank_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<structure_composition<rna4, dot_bracket3>>, uint8_t>));
}

// alphabet_concept: char_type
TEST(structure_composition, char_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<structure_composition<rna4, dot_bracket3>>,
                                underlying_char_t<rna4>>));
}

// alphabet_concept: alphabet_size
TEST(structure_composition, alphabet_size_v)
{
    EXPECT_EQ((alphabet_size_v<structure_composition<rna4, dot_bracket3>>),
              (alphabet_size_v<rna4> * alphabet_size_v<dot_bracket3>));
}

// alphabet_concept: to_rank
TEST(structure_composition, to_rank)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 1);
    EXPECT_EQ(to_rank(std::get<1>(t0)), 2);
    EXPECT_EQ(to_rank(t0),
              to_rank(std::get<0>(t0)) +
              alphabet_size_v<rna4> * to_rank(std::get<1>(t0)));
}

// alphabet_concept: assign_rank
TEST(structure_composition, assign_rank)
{
    using type = structure_composition<rna4, dot_bracket3>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }
}

// alphabet_concept: to_char
TEST(structure_composition, to_char)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), ')');
    EXPECT_EQ(to_char(t0), 'C');
}

// alphabet_concept: assign_char
TEST(structure_composition, assign_char)
{
    using type = structure_composition<rna4, dot_bracket3>;

    type t0{rna4::C, dot_bracket3::PAIR_OPEN};
    char qchar = to_char(std::get<1>(t0));

    assign_char(t0, 'A');
    EXPECT_EQ(to_char(t0), 'A');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'C');
    EXPECT_EQ(to_char(t0), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'G');
    EXPECT_EQ(to_char(t0), 'G');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'U');
    EXPECT_EQ(to_char(t0), 'U');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'N');
    EXPECT_EQ(to_char(t0), 'A');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
}

// alphabet_concept: stream
TEST(structure_composition, outstream)
{
    structure_composition<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};
    std::stringstream s;
    s << t0;
    t0 = rna4::A;
    s << t0;

    EXPECT_EQ(s.str(), "CA");
}
