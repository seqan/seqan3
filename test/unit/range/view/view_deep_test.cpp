// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/take.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

#include "view_concept_check.hpp"

namespace seqan3::view
{
inline auto const deep_reverse = deep{std::view::reverse};
inline auto const deep_take = deep{ranges::view::take};
inline auto const deep_take2 = deep{ranges::view::take(2)};
}

using namespace seqan3;

// ------------------------------------------------------------------
// no parameters
// ------------------------------------------------------------------

TEST(view_deep_reverse, basic)
{
    dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    dna5_vector v0 = foo | view::deep{std::view::reverse};
    EXPECT_EQ(v0, "ATGCA"_dna5);

    // pipe notation
    dna5_vector v = foo | view::deep_reverse;
    EXPECT_EQ(v, "ATGCA"_dna5);

    // function notation
    dna5_vector v2(view::deep_reverse(foo));
    EXPECT_EQ(v2, "ATGCA"_dna5);

    // combinability
    dna5_vector v3 = foo | view::deep_reverse | std::view::reverse;
    EXPECT_EQ(v3, "ACGTA"_dna5);
}

TEST(view_deep_reverse, deep)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | view::deep_reverse;

    ASSERT_EQ(size(v), 2u);
    EXPECT_TRUE((std::ranges::equal(v[0], "ATGCA"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TACGT"_dna5)));
}

TEST(view_deep_reverse, concepts)
{

    using namespace seqan3::test;
    std::vector<dna5_vector> vec{"ACGTA"_dna5, "TGCAT"_dna5};
    auto v1 = vec | view::deep_reverse;
    auto v_elem = v1[0];

    // --- View ---

    EXPECT_TRUE((preserved<decltype(vec), decltype(v1)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v1)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v1)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v1)>({Contiguous, Output})));

    // --- Underlying element ---

    EXPECT_TRUE((preserved<decltype(vec), decltype(v_elem)>({Input, Forward, Bidirectional, RandomAccess, Sized, Common, 
                                                         Output, ConstIterable})));
    EXPECT_TRUE((guaranteed<decltype(vec), decltype(v_elem)>({View})));
    EXPECT_TRUE(weak_guaranteed<decltype(v_elem)>({Viewable}));
    EXPECT_TRUE((lost<decltype(vec), decltype(v_elem)>({Contiguous})));
}

// ------------------------------------------------------------------
// parameters preserved
// ------------------------------------------------------------------

TEST(view_deep_take, basic)
{
    dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    dna5_vector v0 = foo | view::deep{ranges::view::take}(2);
    EXPECT_EQ(v0, "AC"_dna5);

    // pipe notation
    dna5_vector v = foo | view::deep_take(2);
    EXPECT_EQ(v, "AC"_dna5);

    // function notation
    dna5_vector v2(view::deep_take(foo, 2));
    EXPECT_EQ(v2, "AC"_dna5);

    // combinability
    dna5_vector v3 = foo | view::deep_take(2) | std::view::reverse;
    EXPECT_EQ(v3, "CA"_dna5);
}

TEST(view_deep_take, deep)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5, "FOO"_dna5};

    auto v = foo | view::deep_take(2);

    ASSERT_EQ(size(v), 3u);
    EXPECT_TRUE((std::ranges::equal(v[0], "AC"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TG"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[2], "NN"_dna5)));

    int i = 2;
    auto v2 = foo | view::deep_take(i);

    ASSERT_EQ(size(v2), 3u);
    EXPECT_TRUE((std::ranges::equal(v2[0], "AC"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v2[1], "TG"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v2[2], "NN"_dna5)));
}

// ------------------------------------------------------------------
// parameters hardcoded
// ------------------------------------------------------------------

TEST(view_deep_take2, basic)
{
    dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    dna5_vector v0 = foo | view::deep{ranges::view::take(2)};
    EXPECT_EQ(v0, "AC"_dna5);

    // pipe notation
    dna5_vector v = foo | view::deep_take2;
    EXPECT_EQ(v, "AC"_dna5);

    // function notation
    dna5_vector v2(view::deep_take2(foo));
    EXPECT_EQ(v2, "AC"_dna5);

    // combinability
    dna5_vector v3 = foo | view::deep_take2 | std::view::reverse;
    EXPECT_EQ(v3, "CA"_dna5);
}

TEST(view_deep_take2, deep)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5, "FOO"_dna5};

    auto v = foo | view::deep_take2;

    ASSERT_EQ(size(v), 3u);
    EXPECT_TRUE((std::ranges::equal(v[0], "AC"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TG"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[2], "NN"_dna5)));
}
