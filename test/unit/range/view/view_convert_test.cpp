// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;

TEST(view_convert, basic)
{
    std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};
    std::vector<bool> cmp{1, 1, 0, 1, 0, 0, 1, 1, 1};

    // pipe notation
    std::vector<bool> v = vec | view::convert<bool>;
    EXPECT_EQ(cmp, v);

    // function notation
    std::vector<bool> v2(view::convert<bool>(vec));
    EXPECT_EQ(cmp, v2);

    // combinability
    std::vector<bool> cmp2{1, 1, 1, 0, 0, 1, 0, 1, 1};
    std::vector<bool> v3 = vec | view::convert<bool> | std::view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, explicit_conversion)
{
    dna5_vector vec{"ACGNTNGGN"_dna5};
    dna4_vector cmp{"ACGATAGGA"_dna4};

    // pipe notation
    dna4_vector v = vec | view::convert<dna4>;
    EXPECT_EQ(cmp, v);

    // function notation
    dna4_vector v2(view::convert<dna4>(vec));
    EXPECT_EQ(cmp, v2);

    // combinability
    dna4_vector cmp2{"AGGATAGCA"_dna4};
    dna4_vector v3 = vec | view::convert<dna4> | std::view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, concepts)
{
    using namespace seqan3::test;
    dna5_vector vec{"ACGNTNGGN"_dna5};
    auto v1 = vec | view::convert<dna4>;
    
    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v1)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Contiguous, Output})));
}
