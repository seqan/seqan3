// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/trim.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;
using namespace seqan3::view;

TEST(view_trim, standalone)
{
    std::vector<phred42> vec{ phred42{40}, phred42{40}, phred42{30}, phred42{20}, phred42{10}};
    std::vector<phred42> cmp1{phred42{40}, phred42{40}, phred42{30}, phred42{20}};
    std::vector<phred42> cmp2{phred42{40}, phred42{40}};

    // trim by phred_value
    auto v1 = vec | view::trim(20u);                        // == ['I','I','?','5']
    EXPECT_EQ(std::vector<phred42>(v1), cmp1);

    // trim by quality character
    auto v2 = vec | view::trim(phred42{40});             // == ['I','I']
    EXPECT_EQ(std::vector<phred42>(v2), cmp2);

    // function syntax
    auto v3 = view::trim(vec, 20u);                         // == ['I','I','?','5']
    EXPECT_EQ(std::vector<phred42>(v3), cmp1);

    // combinability
    std::string v4 = view::trim(vec, 20u) | view::to_char;  // == "II?5"
    EXPECT_EQ("II?5", v4);
}

TEST(view_trim, qualified)
{
    std::vector<dna5q> vec{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}},
                           {'A'_dna5, phred42{20}}, {'T'_dna5, phred42{10}}};
    std::vector<dna5q> cmp1{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}},
                            {'A'_dna5, phred42{20}}};
    std::vector<dna5q> cmp2{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}};

    // trim by phred_value
    auto v1 = vec | view::trim(20u);
    EXPECT_EQ(std::vector<dna5q>(v1), cmp1);

    // trim by quality character
    auto v2 = vec | view::trim(dna5q{'C'_dna5, phred42{40}});
    EXPECT_EQ(std::vector<dna5q>(v2), cmp2);

    // function syntax
    auto v3 = view::trim(vec, 20u);
    EXPECT_EQ(std::vector<dna5q>(v3), cmp1);

    // combinability
    std::string v4 = view::trim(vec, 20u) | view::to_char;
    EXPECT_EQ("AGGA", v4);
}

TEST(view_trim, concepts)
{
    using namespace seqan3::test;
    std::vector<dna5q> vec{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}}, {'A'_dna5, phred42{20}}, {'T'_dna5, phred42{10}}};
    auto v1 = vec | view::trim(20u);

    //TODO Contiguous should be preserved, but isn't
    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>(
        {Input, Forward, Bidirectional, RandomAccess/*, Contiguous*/, Output, ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v1)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Sized, Common})));
}