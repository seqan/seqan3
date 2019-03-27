// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;

TEST(view_to_rank, basic)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    std::vector<unsigned> cmp{0,1,4,4,4,2,0,4,0};

    // pipe notation
    std::vector<unsigned> v = vec | view::to_rank;
    EXPECT_EQ(cmp, v);

    // function notation
    std::vector<unsigned> v2(view::to_rank(vec));
    EXPECT_EQ(cmp, v2);

    // combinability
    std::vector<unsigned> cmp2{0, 4, 0, 2, 4, 4, 4, 1, 0};
    std::vector<unsigned> v3 = vec | view::to_rank | std::view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_to_rank, concepts)
{
    using namespace seqan3::test;
    dna5_vector vec{"ACTTTGATA"_dna5};
    auto v1 = vec | view::to_rank;

    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v1)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Contiguous, Output})));
}
