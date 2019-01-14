// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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

TEST(view_get, basic)
{
    // TODO remove const-ness from input vector once alphabet_proxy's complement doesnt cause ICE
    std::vector<dna4q> const qv{{'A'_dna4, phred42{0}}, {'C'_dna4, phred42{1}}, {'G'_dna4, phred42{2}}, {'T'_dna4, phred42{3}}};
    std::vector<dna4> cmp0{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
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
    // TODO remove const-ness from input vector once alphabet_proxy inherits it's alphabet
    std::vector<qualified<masked<dna4>, phred42>> const t{{{'A'_dna4, mask::MASKED}, phred42{0}},
                                                    {{'C'_dna4, mask::UNMASKED}, phred42{1}},
                                                    {{'G'_dna4, mask::MASKED}, phred42{2}},
                                                    {{'T'_dna4, mask::UNMASKED}, phred42{3}}};

    // functor notation
    std::vector<masked<dna4>> cmp0{{'A'_dna4, mask::MASKED}, {'C'_dna4, mask::UNMASKED},
                                  {'G'_dna4, mask::MASKED}, {'T'_dna4, mask::UNMASKED}};
    std::vector<masked<dna4>> functor0 = view::get<0>(t);
    EXPECT_EQ(cmp0, functor0);

    std::vector<phred42> cmp1{phred42{0}, phred42{1}, phred42{2}, phred42{3}};
    std::vector<phred42> functor1 = view::get<1>(t);
    EXPECT_EQ(cmp1, functor1);

    std::vector<dna4> cmp00{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
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
    std::vector<masked<dna4>> cmprev{{'T'_dna4, mask::UNMASKED}, {'G'_dna4, mask::MASKED},
                                     {'C'_dna4, mask::UNMASKED}, {'A'_dna4, mask::MASKED}};
    std::vector<masked<dna4>> revtest = t | view::get<0> | view::reverse;
    EXPECT_EQ(cmprev, revtest);

    std::vector<dna4> cmprev2{'T'_dna4, 'G'_dna4, 'C'_dna4, 'A'_dna4};
    std::vector<dna4> revtest2 = t | view::get<0> | view::get<0> | view::reverse;
    EXPECT_EQ(cmprev2, revtest2);

    // reference check
    functor0[0] = masked<dna4>{'T'_dna4, mask::UNMASKED};
    std::vector<masked<dna4>> cmpref{{'T'_dna4, mask::UNMASKED}, {'C'_dna4, mask::UNMASKED},
                                     {'G'_dna4, mask::MASKED}, {'T'_dna4, mask::UNMASKED}};
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
