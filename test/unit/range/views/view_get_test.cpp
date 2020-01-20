// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_get, basic)
{
    // TODO remove const-ness from input vector once alphabet_proxy's complement doesnt cause ICE
    std::vector<dna4q> const qv{{'A'_dna4, phred42{0}}, {'C'_dna4, phred42{1}}, {'G'_dna4, phred42{2}}, {'T'_dna4, phred42{3}}};
    std::vector<dna4> cmp0{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    std::vector<phred42> cmp1{phred42{0}, phred42{1}, phred42{2}, phred42{3}};

    //functor
    dna4_vector functor0 = views::get<0>(qv) | views::to<std::vector>;
    std::vector<phred42> functor1 = views::get<1>(qv) | views::to<std::vector>;
    EXPECT_EQ(cmp0, functor0);
    EXPECT_EQ(cmp1, functor1);

    // pipe notation
    dna4_vector pipe0 = qv | views::get<0> | views::to<std::vector>;
    std::vector<phred42> pipe1 = qv | views::get<1> | views::to<std::vector>;
    EXPECT_EQ(cmp0, pipe0);
    EXPECT_EQ(cmp1, pipe1);

    // combinability
    dna4_vector cmp2{"TGCA"_dna4};
    dna4_vector comp = qv | views::get<0> | views::complement | views::to<std::vector>;
    EXPECT_EQ(cmp2, comp);

    std::string cmp3{"TGCA"};
    std::string to_char_test = comp | views::to_char | views::to<std::string>;
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
    std::vector<masked<dna4>> functor0 = views::get<0>(t) | views::to<std::vector>;
    EXPECT_EQ(cmp0, functor0);

    std::vector<phred42> cmp1{phred42{0}, phred42{1}, phred42{2}, phred42{3}};
    std::vector<phred42> functor1 = views::get<1>(t) | views::to<std::vector>;
    EXPECT_EQ(cmp1, functor1);

    std::vector<dna4> cmp00{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    std::vector<dna4> functor00 = views::get<0>(views::get<0>(t)) | views::to<std::vector>;
    EXPECT_EQ(cmp00, functor00);

    // pipe notation
    std::vector<masked<dna4>> pipe0 = t | views::get<0> | views::to<std::vector>;
    EXPECT_EQ(cmp0, pipe0);

    std::vector<phred42> pipe1 = t | views::get<1> | views::to<std::vector>;
    EXPECT_EQ(cmp1, pipe1);

    std::vector<dna4> pipe00 = t | views::get<0> | views::get<0> | views::to<std::vector>;
    EXPECT_EQ(cmp00, pipe00);

    // combinability
    std::vector<masked<dna4>> cmprev{{'T'_dna4, mask::UNMASKED}, {'G'_dna4, mask::MASKED},
                                     {'C'_dna4, mask::UNMASKED}, {'A'_dna4, mask::MASKED}};
    std::vector<masked<dna4>> revtest = t | views::get<0> | std::views::reverse | views::to<std::vector>;
    EXPECT_EQ(cmprev, revtest);

    std::vector<dna4> cmprev2{'T'_dna4, 'G'_dna4, 'C'_dna4, 'A'_dna4};
    std::vector<dna4> revtest2 = t | views::get<0> | views::get<0> | std::views::reverse | views::to<std::vector>;
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
    std::vector<int> pair_func = views::get<0>(pair_test) | views::to<std::vector>;
    std::vector<int> tuple_func = views::get<0>(tuple_test) | views::to<std::vector>;
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
    std::vector<int> pair_pipe = pair_test | views::get<0> | views::to<std::vector>;
    std::vector<int> tuple_pipe = tuple_test | views::get<0> | views::to<std::vector>;
    EXPECT_EQ(cmp, pair_pipe);
    EXPECT_EQ(cmp, tuple_pipe);
}

TEST(view_get, concepts)
{
    std::vector<std::tuple<int, int>> vec{{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), std::tuple<int, int>>));

    auto v1 = vec | views::get<0>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), std::tuple<int, int>>));
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), int>));
}

// https://github.com/seqan/seqan3/issues/745
TEST(view_get, nested_zip_view)
{
    std::vector vec1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    auto get_view = views::zip(views::zip(vec1, vec1), vec1) | views::get<0>;

    for (auto && elem : get_view)
        std::get<0>(elem) = -1;

    EXPECT_EQ(vec1[0], -1);
}
