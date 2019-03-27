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
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;

TEST(view_to_char, basic)
{
    dna5_vector vec{"ACTTTGATA"_dna5};
    std::string cmp{"ACTTTGATA"};

    // pipe notation
    std::string v = vec | view::to_char;
    EXPECT_EQ(cmp, v);

    // function notation
    std::string v2(view::to_char(vec));
    EXPECT_EQ(cmp, v2);

    // combinability
    std::string cmp2{"ATAGTTTCA"};
    std::string v3 = vec | view::to_char | std::view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_to_char, concepts)
{
    using namespace seqan3::test;
    dna5_vector vec{"ACTTTGATA"_dna5};
    auto v1 = vec | view::to_char;

    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v1)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Contiguous, Output})));
}
