// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <iostream>
#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#ifdef SEQAN3_DEPRECATED_310
#include <seqan3/range/views/move.hpp>
#endif // SEQAN3_DEPRECATED_310
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

using seqan3::operator""_dna5;

TEST(view_move, basic)
{
    std::string vec{"ACTTTGATA"};

    // pipe notation
    auto v = vec | seqan3::views::move;
    EXPECT_RANGE_EQ(vec, v); // equality comparison does not move

    // function notation
    auto v2(seqan3::views::move(vec));
    EXPECT_RANGE_EQ(vec, v2); // equality comparison does not move

    // combinability
    seqan3::dna5_vector vec2{"ACGTA"_dna5};

    EXPECT_RANGE_EQ("TGCAT"_dna5, vec2 | seqan3::views::complement
                                       | seqan3::views::move /*NOP, because already temporaries*/);
}

TEST(view_move, concepts)
{
    seqan3::dna5_vector vec{"ACTTTGATA"_dna5};
    auto v1 = vec | seqan3::views::move;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), char>));

    EXPECT_TRUE((std::is_same_v<decltype(v1[0]), seqan3::dna5 &&>));

    // complement generates values
    auto v2 = vec | seqan3::views::complement | seqan3::views::move;
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v2)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), char>));

    EXPECT_TRUE((std::is_same_v<decltype(v2[0]), seqan3::dna5>)); // don't add const-ness to values
}

#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
