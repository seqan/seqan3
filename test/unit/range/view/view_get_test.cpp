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

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/std/view/reverse.hpp>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/mask/masked.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(view_get, basic)
{
    std::vector<dna4q> qv{{dna4::A, phred42{0}}, {dna4::C, phred42{1}}, {dna4::G, phred42{2}}, {dna4::T, phred42{3}}};
    std::vector<dna4> cmp0{dna4::A, dna4::C, dna4::G, dna4::T};
    std::vector<phred42> cmp1{phred42{0}, phred42{1}, phred42{2}, phred42{3}};

    //functor
    dna4_vector functor0 = view::get<0>(qv);
    std::vector<phred42> functor1 = view::get<1>(qv);
    EXPECT_EQ(cmp0, functor0);
    EXPECT_EQ(cmp1, functor1);

    // pipe notation
    dna4_vector pipe0 = qv | view::get<0>;
    std::vector<phred42> pipe1 = qv | view::get<1>;
    EXPECT_EQ(cmp0, pipe0);
    EXPECT_EQ(cmp1, pipe1);

    // combinability
    dna4_vector cmp2{"TGCA"_dna4};
    dna4_vector comp = qv | view::get<0> | view::complement;
    EXPECT_EQ(cmp2, comp);

    std::string cmp3{"TGCA"};
    std::string to_char_test = comp | view::to_char;
    EXPECT_EQ(cmp3, to_char_test);

    // reference return check
    functor1[0] = phred42{4};
    std::vector<phred42> cmp4{phred42{4}, phred42{1}, phred42{2}, phred42{3}};
    EXPECT_EQ(cmp4, functor1);
}

TEST(view_get, advanced)
{
    std::vector<qualified<masked<dna4>, phred42>> t{{{dna4::A, mask::MASKED}, phred42{0}},
                                                    {{dna4::C, mask::UNMASKED}, phred42{1}},
                                                    {{dna4::G, mask::MASKED}, phred42{2}},
                                                    {{dna4::T, mask::UNMASKED}, phred42{3}}};

    // functor notation
    std::vector<masked<dna4>> cmp0{{dna4::A, mask::MASKED}, {dna4::C, mask::UNMASKED},
                                  {dna4::G, mask::MASKED}, {dna4::T, mask::UNMASKED}};
    std::vector<masked<dna4>> functor0 = view::get<0>(t);
    EXPECT_EQ(cmp0, functor0);

    std::vector<phred42> cmp1{phred42{0}, phred42{1}, phred42{2}, phred42{3}};
    std::vector<phred42> functor1 = view::get<1>(t);
    EXPECT_EQ(cmp1, functor1);

    std::vector<dna4> cmp00{dna4::A, dna4::C, dna4::G, dna4::T};
    std::vector<dna4> functor00 = view::get<0>(view::get<0>(t));
    EXPECT_EQ(cmp00, functor00);

    // pipe notation
    std::vector<masked<dna4>> pipe0 = t | view::get<0>;
    EXPECT_EQ(cmp0, pipe0);

    std::vector<phred42> pipe1 = t | view::get<1>;
    EXPECT_EQ(cmp1, pipe1);

    std::vector<dna4> pipe00 = t | view::get<0> | view::get<0>;
    EXPECT_EQ(cmp00, pipe00);

    // combinability
    std::vector<masked<dna4>> cmprev{{dna4::T, mask::UNMASKED}, {dna4::G, mask::MASKED},
                                     {dna4::C, mask::UNMASKED}, {dna4::A, mask::MASKED}};
    std::vector<masked<dna4>> revtest = t | view::get<0> | view::reverse;
    EXPECT_EQ(cmprev, revtest);

    std::vector<dna4> cmprev2{dna4::T, dna4::G, dna4::C, dna4::A};
    std::vector<dna4> revtest2 = t | view::get<0> | view::get<0> | view::reverse;
    EXPECT_EQ(cmprev2, revtest2);

    // reference check
    functor0[0] = masked<dna4>{dna4::T, mask::UNMASKED};
    std::vector<masked<dna4>> cmpref{{dna4::T, mask::UNMASKED}, {dna4::C, mask::UNMASKED},
                                     {dna4::G, mask::MASKED}, {dna4::T, mask::UNMASKED}};
    EXPECT_EQ(cmpref, functor0);
}

TEST(view_get, tuple_pair)
{
    std::vector<std::pair<int, int>> pair_test{{0, 1}, {1, 2}, {2, 3}, {3, 4}};
    std::vector<std::tuple<int, int>> tuple_test{{0, 1}, {1, 2}, {2, 3}, {3, 4}};

    // functor notation
    std::vector<int> cmp{0, 1, 2, 3};
    std::vector<int> pair_func = view::get<0>(pair_test);
    std::vector<int> tuple_func = view::get<0>(tuple_test);
    EXPECT_EQ(cmp, pair_func);
    EXPECT_EQ(cmp, tuple_func);

    // reference test
    cmp[0] = 4;
    pair_func[0] = 4;
    tuple_func[0] = 4;
    EXPECT_EQ(cmp, pair_func);
    EXPECT_EQ(cmp, tuple_func);

    // pipe notation
    cmp[0] = 0;
    std::vector<int> pair_pipe = pair_test | view::get<0>;
    std::vector<int> tuple_pipe = tuple_test | view::get<0>;
    EXPECT_EQ(cmp, pair_pipe);
    EXPECT_EQ(cmp, tuple_pipe);
}

TEST(view_get, concepts)
{
    std::vector<std::tuple<int, int>> vec{{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), std::tuple<int, int>>));

    auto v1 = vec | view::get<0>;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), std::tuple<int, int>>));
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), int>));
}
