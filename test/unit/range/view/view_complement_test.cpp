// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;

TEST(view_complement, basic)
{
    dna5_vector foo{"ACGTA"_dna5};

    // pipe notation
    dna5_vector v = foo | view::complement;
    EXPECT_EQ(v, "TGCAT"_dna5);

    // function notation
    dna5_vector v2(view::complement(foo));
    EXPECT_EQ(v2, "TGCAT"_dna5);

    // combinability
    dna5_vector v3 = foo | view::complement | std::view::reverse;
    EXPECT_EQ(v3, "TACGT"_dna5);
}

TEST(view_complement, deep_view)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | view::complement;

    ASSERT_EQ(size(v), 2u);
    EXPECT_TRUE((std::ranges::equal(v[0], "TGCAT"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "ACGTA"_dna5)));
}

TEST(view_complement, concepts)
{
    using namespace seqan3::test;
    dna5_vector vec{"ACGTA"_dna5};
    auto v1 = vec | view::complement;

    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v1)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Contiguous, Output})));
}
