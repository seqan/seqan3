// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

using namespace seqan3;

/************** ALPHABET and QUALITY concept **********************/

TEST(qualified, concept)
{
    // if the first template argument is a nucleotide, qualified models the nucleotide concept
    EXPECT_TRUE((nucleotide_concept<qualified<dna4, phred42>>));
    // else, qualified models the nucleotide concept
    EXPECT_TRUE((alphabet_concept<qualified<char, phred42>>));
}

TEST(qualified, rank_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<qualified<dna4, phred42>>,
                               uint8_t>));
}

TEST(qualified, char_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<qualified<dna4, phred42>>,
                               underlying_char_t<dna4>>));
}

TEST(qualified, phred_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_phred_t<qualified<dna4, phred42>>,
                               underlying_phred_t<phred42>>));
}

TEST(qualified, alphabet_size_v)
{
    EXPECT_EQ((alphabet_size_v<qualified<dna4, phred42>>),
              (alphabet_size_v<dna4> * alphabet_size_v<phred42>));
}

TEST(qualified, to_rank)
{
    qualified<dna4, phred42> t0{dna4::C, phred42{6}};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 1);
    EXPECT_EQ(to_rank(std::get<1>(t0)), 6);
    EXPECT_EQ(to_rank(t0),
              to_rank(std::get<0>(t0)) +
              alphabet_size_v<dna4> * to_rank(std::get<1>(t0)));
}

TEST(qualified, assign_rank)
{
    using type = qualified<dna4, phred42>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }
}

TEST(qualified, to_char)
{
    qualified<dna4, phred42> t0{dna4::C, phred42{6}};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), '!' + 6);
    EXPECT_EQ(to_char(t0), 'C');
}

TEST(qualified, assign_char)
{
    using type = qualified<dna4, phred42>;

    type t0{dna4::C, phred42{17}};
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
    assign_char(t0, 'T');
    EXPECT_EQ(to_char(t0), 'T');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'N');
    EXPECT_EQ(to_char(t0), 'A');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
}

TEST(qualified, to_phred)
{
    qualified<dna4, phred42> t0{dna4::C, phred42{6}};
    EXPECT_EQ(to_phred(std::get<1>(t0)), 6);
    EXPECT_EQ(to_phred(t0), 6);
}

TEST(qualified, assign_phred)
{
    using type = qualified<dna4, phred42>;

    type t0{dna4::C, phred42{17}};
    char schar = to_char(t0);

    assign_phred(t0, 12);
    EXPECT_EQ(to_phred(t0), 12);
    EXPECT_EQ(to_char(t0), schar);
    assign_phred(t0, 37);
    EXPECT_EQ(to_phred(t0), 37);
    EXPECT_EQ(to_char(t0), schar);
}

TEST(qualified, outstream)
{
    qualified<dna4, phred42> t0{dna4::C, phred42{6}};
    std::stringstream s;
    s << t0;
    t0 = dna4::A;
    s << t0;

    EXPECT_EQ(s.str(), "CA");
}

TEST(qualified, complement)
{
    qualified<dna4, phred42> t0{dna4::A, phred42{8}};
    qualified<dna4, phred42> t0_c{(dna4::A).complement(), phred42{8}};

    EXPECT_EQ(t0.complement(), t0_c);
}
