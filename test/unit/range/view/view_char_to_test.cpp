// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

using namespace seqan3;

TEST(view_char_to, basic)
{
    std::string vec{"ACTTTGATA"};
    dna5_vector cmp{"ACTTTGATA"_dna5};

    // pipe notation
    dna5_vector v = vec | view::char_to<dna5>;
    EXPECT_EQ(cmp, v);

    // function notation
    dna5_vector v2(view::char_to<dna5>(vec));
    EXPECT_EQ(cmp, v2);

    // combinability
    dna5_vector cmp2{"ATAGTTTCA"_dna5};
    dna5_vector v3 = vec | view::char_to<dna5> | std::view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_char_to, deep_view)
{
    std::vector<std::string> foo{"ACGTA", "TGCAT"};

    std::vector<dna5_vector> v = foo | view::char_to<dna5>;

    ASSERT_EQ(size(v), 2u);
    EXPECT_TRUE((std::ranges::equal(v[0], "ACGTA"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TGCAT"_dna5)));
}

TEST(view_char_to, concepts)
{
    using namespace seqan3::test;
    
    std::string vec{"ACTTTGATA"};
    auto v1 = vec | view::char_to<dna5>;

    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View, Viewable})));
    EXPECT_TRUE((weak_guaranteed<decltype(v1)>({Viewable})));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Contiguous, Output})));
}
