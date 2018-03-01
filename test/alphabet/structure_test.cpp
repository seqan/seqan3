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
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

template <typename T>
class structure : public ::testing::Test
{};

// add all alphabets from the structure sub module here
using structure_types  = ::testing::Types<dot_bracket3, wuss51, dssp9>;

TYPED_TEST_CASE(structure, structure_types);

// assign_char functions
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
    else if constexpr (std::is_same_v<TypeParam, wuss51>)
    {
        cmp =
        {
            t::UNPAIRED, t::PAIR_OPEN1, t::PAIR_CLOSE1,
            t::UNPAIRED1, t::UNPAIRED2, t::UNPAIRED3, t::UNPAIRED4, t::UNPAIRED5, t::UNPAIRED6,
            t::PAIR_OPEN, t::PAIR_CLOSE, t::PAIR_OPEN2, t::PAIR_CLOSE2, t::PAIR_OPEN3, t::PAIR_CLOSE3,
            "H"_wuss51.front(),
            "B"_wuss51.front(),
            "E"_wuss51.front(),
            "G"_wuss51.front(),
            "I"_wuss51.front(),
            "T"_wuss51.front(),
            "S"_wuss51.front()
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

// to_char functions
TYPED_TEST(structure, to_char)
{
    if constexpr (std::is_same_v<TypeParam, dot_bracket3>)
    {
        EXPECT_EQ(to_char(TypeParam::UNPAIRED), '.');
        EXPECT_EQ(to_char(TypeParam::PAIR_OPEN), '(');
        EXPECT_EQ(to_char(TypeParam::PAIR_CLOSE), ')');
        EXPECT_EQ(to_char(TypeParam::UNKNOWN), '.');
    }
    else if constexpr (std::is_same_v<TypeParam, wuss51>)
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

// ------------------------------------------------------------------
// concepts
// ------------------------------------------------------------------

TYPED_TEST(structure, concept)
{
    EXPECT_TRUE(alphabet_concept<TypeParam>);
}

TYPED_TEST(structure, rna_structure_concept)
{
    EXPECT_TRUE(rna_structure_concept<dot_bracket3>);
    EXPECT_TRUE(rna_structure_concept<wuss51>);
    EXPECT_TRUE(rna_structure_concept<wuss<>>);  // also wuss51
    EXPECT_TRUE(rna_structure_concept<wuss<67>>);
}

// ------------------------------------------------------------------
// stream operator
// ------------------------------------------------------------------

TEST(structure_stream_operator, dot_bracket3)
{
    std::stringstream ss;
    ss << dot_bracket3::PAIR_OPEN << dot_bracket3::PAIR_CLOSE << dot_bracket3::UNPAIRED;
    EXPECT_EQ(ss.str(), "().");
}

TEST(structure_stream_operator, wuss)
{
    std::stringstream ss;
    ss << wuss51::PAIR_OPEN << wuss51::PAIR_CLOSE << wuss51::UNPAIRED;
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
    EXPECT_EQ(vec1, "((((("_db3);

    std::vector<dot_bracket3> vec2{dot_bracket3::UNPAIRED, dot_bracket3::PAIR_OPEN, dot_bracket3::PAIR_OPEN,
                                   dot_bracket3::PAIR_CLOSE, dot_bracket3::PAIR_CLOSE, dot_bracket3::UNKNOWN};
    EXPECT_EQ(vec2, ".(())."_db3);
}

TEST(wuss_literals, vector)
{
    std::vector<wuss51> vec1;
    vec1.resize(5, wuss51::PAIR_OPEN);
    EXPECT_EQ(vec1, "<<<<<"_wuss51);

    std::vector<wuss51> vec2{wuss51::UNPAIRED, wuss51::PAIR_OPEN, wuss51::PAIR_OPEN,
                             wuss51::PAIR_CLOSE, wuss51::PAIR_CLOSE, wuss51::UNPAIRED};
    EXPECT_EQ(vec2, ".<<>>."_wuss51);
}

TEST(dssp9_literals, vector)
{
    std::vector<dssp9> vec1;
    vec1.resize(5, dssp9::H);
    EXPECT_EQ(vec1, "HHHHH"_dssp9);

    std::vector<dssp9> vec2{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
    EXPECT_EQ(vec2, "EHHHTG"_dssp9);
}

// ------------------------------------------------------------------
// rna_structure
// ------------------------------------------------------------------

TEST(rna_structure, dot_bracket3)
{
    EXPECT_FALSE(dot_bracket3::pseudoknot_support);
    EXPECT_TRUE(dot_bracket3::UNPAIRED.is_unpaired());
    EXPECT_TRUE(dot_bracket3::PAIR_OPEN.is_pair_open());
    EXPECT_TRUE(dot_bracket3::PAIR_CLOSE.is_pair_close());
    EXPECT_TRUE(dot_bracket3::UNKNOWN.is_unpaired());
}

TEST(rna_structure, wuss)
{
    EXPECT_TRUE(wuss51::pseudoknot_support);
    std::vector<wuss51> vec = ".:,-_~;<>()[]{}AaBbCcDd"_wuss51;
    for (unsigned idx = 0; idx <= 6; ++idx)
    {
        EXPECT_TRUE(vec[idx].is_unpaired());
        EXPECT_FALSE(vec[idx].is_pair_open());
        EXPECT_FALSE(vec[idx].is_pair_close());
    }
    for (unsigned idx = 7; idx <= 21; idx+=2)
    {
        EXPECT_TRUE(vec[idx].is_pair_open());
        EXPECT_FALSE(vec[idx].is_unpaired());
        EXPECT_FALSE(vec[idx].is_pair_close());
    }
    for (unsigned idx = 8; idx <= 22; idx+=2)
    {
        EXPECT_TRUE(vec[idx].is_pair_close());
        EXPECT_FALSE(vec[idx].is_unpaired());
        EXPECT_FALSE(vec[idx].is_pair_open());
    }
}

// ------------------------------------------------------------------
// composition nucleotide x structure
// ------------------------------------------------------------------

// default/zero construction
TEST(structured_rna, ctr)
{
    [[maybe_unused]] structured_rna<rna4, dot_bracket3> t1;
}

// aggregate initialization
TEST(structured_rna, aggr)
{
    [[maybe_unused]]  structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_CLOSE};
}

// zero initialization
TEST(structured_rna, zro)
{
    structured_rna<rna4, dot_bracket3> t1{rna4::A, dot_bracket3::UNKNOWN};
    structured_rna<rna4, dot_bracket3> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(structured_rna, cp_ctr)
{
    structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_OPEN};
    structured_rna<rna4, dot_bracket3> t2{t1};
    structured_rna<rna4, dot_bracket3> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(structured_rna, mv_ctr)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::UNPAIRED};
    structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::UNPAIRED};
    structured_rna<rna4, dot_bracket3> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    structured_rna<rna4, dot_bracket3> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(structured_rna, cp_assgn)
{
    structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::UNPAIRED};
    structured_rna<rna4, dot_bracket3> t2;
    structured_rna<rna4, dot_bracket3> t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(structured_rna, mv_assgn)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::UNPAIRED};
    structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::UNPAIRED};
    structured_rna<rna4, dot_bracket3> t2;
    structured_rna<rna4, dot_bracket3> t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(structured_rna, swap)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};
    structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_OPEN};
    structured_rna<rna4, dot_bracket3> t2{};
    structured_rna<rna4, dot_bracket3> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// seqan3::get<1>
TEST(structured_rna, get_i)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), rna4 &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), dot_bracket3 &>);

    EXPECT_EQ(seqan3::get<0>(t0), rna4::C);
    EXPECT_EQ(seqan3::get<1>(t0), dot_bracket3{dot_bracket3::PAIR_OPEN});
}

// std::get<1>
TEST(structured_rna, stdget_i)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::UNPAIRED};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), rna4 &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), dot_bracket3 &>);

    EXPECT_EQ(std::get<0>(t0), rna4::C);
    EXPECT_EQ(std::get<1>(t0), dot_bracket3{dot_bracket3::UNPAIRED});
}

// structured_rna bindings
TEST(structured_rna, struct_binding)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};
    auto [ i, l ] = t0;

    static_assert(std::is_same_v<decltype(i), rna4>);
    static_assert(std::is_same_v<decltype(l), dot_bracket3>);

    EXPECT_EQ(i, rna4::C);
    EXPECT_EQ(l, dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// seqan3::get<type>
TEST(structured_rna, get_type)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};

    EXPECT_EQ(seqan3::get<rna4>(t0), rna4::C);
    EXPECT_EQ(seqan3::get<dot_bracket3>(t0), dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// std::get<type>
TEST(structured_rna, stdget_type)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};

    EXPECT_EQ(std::get<rna4>(t0), rna4::C);
    EXPECT_EQ(std::get<dot_bracket3>(t0), dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// std::tuple_element
TEST(structured_rna, tuple_element)
{
    using pt = structured_rna<rna4, dot_bracket3>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, rna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, dot_bracket3>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// type deduction
TEST(structured_rna, type_deduce)
{
    structured_rna t0{rna4::C, dot_bracket3{dot_bracket3::PAIR_CLOSE}};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, rna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, dot_bracket3>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// explicit cast to element
TEST(structured_rna, cast_to_element)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};

    auto d = static_cast<rna4>(t0);
    auto q = static_cast<dot_bracket3>(t0);
    static_assert(std::is_same_v<decltype(d), rna4>);
    static_assert(std::is_same_v<decltype(q), dot_bracket3>);

    EXPECT_EQ(d, rna4::C);
    EXPECT_EQ(q, dot_bracket3{dot_bracket3::PAIR_CLOSE});
}

// comparison operators
TEST(structured_rna, cmp)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};
    structured_rna<rna4, dot_bracket3> t1{rna4::C, dot_bracket3::PAIR_CLOSE};
    structured_rna<rna4, dot_bracket3> t2{rna4::G, dot_bracket3::PAIR_CLOSE};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

// alphabet_concept: rank_type
TEST(structured_rna, rank_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<structured_rna<rna4, dot_bracket3>>, uint8_t>));
}

// alphabet_concept: char_type
TEST(structured_rna, char_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<structured_rna<rna4, dot_bracket3>>,
                                underlying_char_t<rna4>>));
}

// alphabet_concept: alphabet_size
TEST(structured_rna, alphabet_size_v)
{
    EXPECT_EQ((alphabet_size_v<structured_rna<rna4, dot_bracket3>>),
              (alphabet_size_v<rna4> * alphabet_size_v<dot_bracket3>));
}

// alphabet_concept: to_rank
TEST(structured_rna, to_rank)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 1);
    EXPECT_EQ(to_rank(std::get<1>(t0)), 2);
    EXPECT_EQ(to_rank(t0),
              to_rank(std::get<0>(t0)) +
              alphabet_size_v<rna4> * to_rank(std::get<1>(t0)));
}

// alphabet_concept: assign_rank
TEST(structured_rna, assign_rank)
{
    using type = structured_rna<rna4, dot_bracket3>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }
}

// alphabet_concept: to_char
TEST(structured_rna, to_char)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_CLOSE};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), ')');
    EXPECT_EQ(to_char(t0), 'C');
}

// alphabet_concept: assign_char
TEST(structured_rna, assign_char)
{
    using type = structured_rna<rna4, dot_bracket3>;

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
TEST(structured_rna, outstream)
{
    structured_rna<rna4, dot_bracket3> t0{rna4::C, dot_bracket3::PAIR_OPEN};
    std::stringstream s;
    s << t0;
    t0 = rna4::A;
    s << t0;

    EXPECT_EQ(s.str(), "CA");
}

// nucleotide concept: complement
TEST(structured_rna, complement)
{
    using type = structured_rna<rna5, dot_bracket3>;
    type tU{rna5::U, dot_bracket3::PAIR_OPEN};
    type tG{rna5::G, dot_bracket3::PAIR_OPEN};
    type tC{rna5::C, dot_bracket3::PAIR_OPEN};
    type tA{rna5::A, dot_bracket3::PAIR_OPEN};
    type tN{rna5::N, dot_bracket3::PAIR_OPEN};
    type tt{rna5::A, dot_bracket3::PAIR_CLOSE};

    // check whether characters are converted correctly
    EXPECT_EQ(to_char(complement(tU)), 'A');
    EXPECT_EQ(to_char(complement(tG)), 'C');
    EXPECT_EQ(to_char(complement(tC)), 'G');
    EXPECT_EQ(to_char(complement(tA)), 'U');
    EXPECT_EQ(to_char(complement(tN)), 'N');

    // check whether structure character is not modified
    EXPECT_EQ(complement(tU), tA);
    EXPECT_EQ(complement(tG), tC);
    EXPECT_EQ(complement(tN), tN);
    EXPECT_NE(complement(tU), tt);

    // complement combinations
    EXPECT_EQ(complement(complement(tU)), tU);
    EXPECT_EQ(complement(tU), complement(tU));
}

// ------------------------------------------------------------------
// composition aminoacid x protein structure
// ------------------------------------------------------------------

// default/zero construction
TEST(structured_aa, ctr)
{
    [[maybe_unused]] structured_aa<aa27, dssp9> t1;
}

// aggregate initialization
TEST(structured_aa, aggr)
{
    [[maybe_unused]] structured_aa<aa27, dssp9> t1{aa27::C, dssp9::B};
}

// zero initialization
TEST(structured_aa, zro)
{
    structured_aa<aa27, dssp9> t1{aa27::A, dssp9::H};
    structured_aa<aa27, dssp9> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(structured_aa, cp_ctr)
{
    structured_aa<aa27, dssp9> t1{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t2{t1};
    structured_aa<aa27, dssp9> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(structured_aa, mv_ctr)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t1{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    structured_aa<aa27, dssp9> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(structured_aa, cp_assgn)
{
    structured_aa<aa27, dssp9> t1{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t2;
    structured_aa<aa27, dssp9> t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(structured_aa, mv_assgn)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t1{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t2;
    structured_aa<aa27, dssp9> t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(structured_aa, swap)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::E};
    structured_aa<aa27, dssp9> t1{aa27::C, dssp9::E};
    structured_aa<aa27, dssp9> t2{};
    structured_aa<aa27, dssp9> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// seqan3::get<1>
TEST(structured_aa, get_i)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::E};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), aa27 &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), dssp9 &>);

    EXPECT_EQ(seqan3::get<0>(t0), aa27::C);
    EXPECT_EQ(seqan3::get<1>(t0), dssp9::E);
}

// std::get<1>
TEST(structured_aa, stdget_i)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), aa27 &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), dssp9 &>);

    EXPECT_EQ(std::get<0>(t0), aa27::C);
    EXPECT_EQ(std::get<1>(t0), dssp9::G);
}

// structured_aa bindings
TEST(structured_aa, struct_binding)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};
    auto [ i, l ] = t0;

    static_assert(std::is_same_v<decltype(i), aa27>);
    static_assert(std::is_same_v<decltype(l), dssp9>);

    EXPECT_EQ(i, aa27::C);
    EXPECT_EQ(l, dssp9::G);
}

// seqan3::get<type>
TEST(structured_aa, get_type)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};

    EXPECT_EQ(seqan3::get<aa27>(t0), aa27::C);
    EXPECT_EQ(seqan3::get<dssp9>(t0), dssp9::G);
}

// std::get<type>
TEST(structured_aa, stdget_type)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};

    EXPECT_EQ(std::get<aa27>(t0), aa27::C);
    EXPECT_EQ(std::get<dssp9>(t0), dssp9::G);
}

// std::tuple_element
TEST(structured_aa, tuple_element)
{
    using pt = structured_aa<aa27, dssp9>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, aa27>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, dssp9>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// type deduction
TEST(structured_aa, type_deduce)
{
    structured_aa t0{aa27::C, dssp9::G};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, aa27>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, dssp9>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// explicit cast to element
TEST(structured_aa, cast_to_element)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};

    auto d = static_cast<aa27>(t0);
    auto q = static_cast<dssp9>(t0);
    static_assert(std::is_same_v<decltype(d), aa27>);
    static_assert(std::is_same_v<decltype(q), dssp9>);

    EXPECT_EQ(d, aa27::C);
    EXPECT_EQ(q, dssp9::G);
}

// comparison operators
TEST(structured_aa, cmp)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::H};
    structured_aa<aa27, dssp9> t1{aa27::C, dssp9::B};
    structured_aa<aa27, dssp9> t2{aa27::G, dssp9::E};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

// alphabet_concept: rank_type
TEST(structured_aa, rank_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<structured_aa<aa27, dssp9>>, uint8_t>));
}

// alphabet_concept: char_type
TEST(structured_aa, char_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<structured_aa<aa27, dssp9>>,
        underlying_char_t<aa27>>));
}

// alphabet_concept: alphabet_size
TEST(structured_aa, alphabet_size_v)
{
    EXPECT_EQ((alphabet_size_v<structured_aa<aa27, dssp9>>),
              (alphabet_size_v<aa27> * alphabet_size_v<dssp9>));
}

// alphabet_concept: to_rank
TEST(structured_aa, to_rank)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 2);
    EXPECT_EQ(to_rank(std::get<1>(t0)), 3);
    EXPECT_EQ(to_rank(t0),
              to_rank(std::get<0>(t0)) +
              alphabet_size_v<aa27> * to_rank(std::get<1>(t0)));
}

// alphabet_concept: assign_rank
TEST(structured_aa, assign_rank)
{
    using type = structured_aa<aa27, dssp9>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }
}

// alphabet_concept: to_char
TEST(structured_aa, to_char)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::G};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), 'G');
    EXPECT_EQ(to_char(t0), 'C');
}

// alphabet_concept: assign_char
TEST(structured_aa, assign_char)
{
    using type = structured_aa<aa27, dssp9>;

    type t0{aa27::C, dssp9::T};
    char qchar = to_char(std::get<1>(t0));

    assign_char(t0, 'A');
    EXPECT_EQ(to_char(t0), 'A');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'Y');
    EXPECT_EQ(to_char(t0), 'Y');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'W');
    EXPECT_EQ(to_char(t0), 'W');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'D');
    EXPECT_EQ(to_char(t0), 'D');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'X');
    EXPECT_EQ(to_char(t0), 'X');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
}

// alphabet_concept: stream
TEST(structured_aa, outstream)
{
    structured_aa<aa27, dssp9> t0{aa27::C, dssp9::T};
    std::stringstream s;
    s << t0;
    t0 = aa27::A;
    s << t0;

    EXPECT_EQ(s.str(), "CA");
}
