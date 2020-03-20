// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/trim.hpp>
#include <seqan3/std/ranges>

using seqan3::operator""_dna5;

TEST(view_trim, standalone)
{
    std::vector<seqan3::phred42> vec{ seqan3::phred42{40}, seqan3::phred42{40},
                                      seqan3::phred42{30}, seqan3::phred42{20}, seqan3::phred42{10}};
    std::vector<seqan3::phred42> cmp1{seqan3::phred42{40}, seqan3::phred42{40},
                                      seqan3::phred42{30}, seqan3::phred42{20}};
    std::vector<seqan3::phred42> cmp2{seqan3::phred42{40}, seqan3::phred42{40}};

    // trim by phred_value
    auto v1 = vec | seqan3::views::trim(20u);                   // == ['I','I','?','5']
    EXPECT_EQ(v1 | seqan3::views::to<std::vector>, cmp1);

    // trim by quality character
    auto v2 = vec | seqan3::views::trim(seqan3::phred42{40});   // == ['I','I']
    EXPECT_EQ(v2 | seqan3::views::to<std::vector>, cmp2);

    // function syntax
    auto v3 = seqan3::views::trim(vec, 20u);                    // == ['I','I','?','5']
    EXPECT_EQ(v3 | seqan3::views::to<std::vector>, cmp1);

    // combinability
    std::string v4 = seqan3::views::trim(vec, 20u) | seqan3::views::to_char | seqan3::views::to<std::string>; //=="II?5"
    EXPECT_EQ("II?5", v4);
}

TEST(view_trim, qualified)
{
    std::vector<seqan3::dna5q> vec{{'A'_dna5, seqan3::phred42{40}},
                                   {'G'_dna5, seqan3::phred42{40}},
                                   {'G'_dna5, seqan3::phred42{30}},
                                   {'A'_dna5, seqan3::phred42{20}},
                                   {'T'_dna5, seqan3::phred42{10}}};
    std::vector<seqan3::dna5q> cmp1{{'A'_dna5, seqan3::phred42{40}},
                                    {'G'_dna5, seqan3::phred42{40}},
                                    {'G'_dna5, seqan3::phred42{30}},
                                    {'A'_dna5, seqan3::phred42{20}}};
    std::vector<seqan3::dna5q> cmp2{{'A'_dna5, seqan3::phred42{40}},
                                    {'G'_dna5, seqan3::phred42{40}}};

    // trim by phred_value
    auto v1 = vec | seqan3::views::trim(20u);
    EXPECT_EQ(v1 | seqan3::views::to<std::vector>, cmp1);

    // trim by quality character
    auto v2 = vec | seqan3::views::trim(seqan3::dna5q{'C'_dna5, seqan3::phred42{40}});
    EXPECT_EQ(v2 | seqan3::views::to<std::vector>, cmp2);

    // function syntax
    auto v3 = seqan3::views::trim(vec, 20u);
    EXPECT_EQ(v3 | seqan3::views::to<std::vector>, cmp1);

    // combinability
    std::string v4 = seqan3::views::trim(vec, 20u) | seqan3::views::to_char | seqan3::views::to<std::string>;
    EXPECT_EQ("AGGA", v4);
}

TEST(view_trim, concepts)
{
    std::vector<seqan3::dna5q> vec{{'A'_dna5, seqan3::phred42{40}},
                                   {'G'_dna5, seqan3::phred42{40}},
                                   {'G'_dna5, seqan3::phred42{30}},
                                   {'A'_dna5, seqan3::phred42{20}},
                                   {'T'_dna5, seqan3::phred42{10}}};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), seqan3::dna5q>));
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);

    auto v1 = vec | seqan3::views::trim(20u);
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), seqan3::dna5q>));
    EXPECT_TRUE(!std::ranges::sized_range<decltype(v1)>);
}
