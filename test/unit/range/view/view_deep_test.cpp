// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <iostream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/take.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/view/reverse.hpp>

namespace seqan3::view
{
inline auto const deep_reverse = deep{view::reverse};
inline auto const deep_take = deep{ranges::view::take};
inline auto const deep_take2 = deep{ranges::view::take(2)};
}

using namespace seqan3;
using namespace seqan3::literal;

// ------------------------------------------------------------------
// no parameters
// ------------------------------------------------------------------

TEST(view_deep_reverse, basic)
{
    dna5_vector foo{"ACGTA"_dna5};

    // pipe notation, temporary
    dna5_vector v0 = foo | view::deep{view::reverse};
    EXPECT_EQ(v0, "ATGCA"_dna5);

    // pipe notation
    dna5_vector v = foo | view::deep_reverse;
    EXPECT_EQ(v, "ATGCA"_dna5);

    // function notation
    dna5_vector v2(view::deep_reverse(foo));
    EXPECT_EQ(v2, "ATGCA"_dna5);

    // combinability
    dna5_vector v3 = foo | view::deep_reverse | view::reverse;
    EXPECT_EQ(v3, "ACGTA"_dna5);
}

TEST(view_deep_reverse, deep)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | view::deep_reverse;

    ASSERT_EQ(size(v), 2);
    EXPECT_TRUE((std::ranges::equal(v[0], "ATGCA"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TACGT"_dna5)));
}

TEST(view_deep_reverse, concepts)
{
    std::vector<dna5_vector> vec{"ACGTA"_dna5, "TGCAT"_dna5};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), dna5_vector>));

    auto v1 = vec | view::deep_reverse;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna5_vector>)); // view temporary returned in deep case

    auto v_elem = v1[0];
    EXPECT_TRUE(std::ranges::InputRange<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::View<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v_elem)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v_elem)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v_elem)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v_elem), dna5>));
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
    dna5_vector v3 = foo | view::deep_take(2) | view::reverse;
    EXPECT_EQ(v3, "CA"_dna5);
}

TEST(view_deep_take, deep)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5, "FOO"_dna5};

    auto v = foo | view::deep_take(2);

    ASSERT_EQ(size(v), 3);
    EXPECT_TRUE((std::ranges::equal(v[0], "AC"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TG"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[2], "NN"_dna5)));

    int i = 2;
    auto v2 = foo | view::deep_take(i);

    ASSERT_EQ(size(v2), 3);
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
    dna5_vector v3 = foo | view::deep_take2 | view::reverse;
    EXPECT_EQ(v3, "CA"_dna5);
}

TEST(view_deep_take2, deep)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5, "FOO"_dna5};

    auto v = foo | view::deep_take2;

    ASSERT_EQ(size(v), 3);
    EXPECT_TRUE((std::ranges::equal(v[0], "AC"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[1], "TG"_dna5)));
    EXPECT_TRUE((std::ranges::equal(v[2], "NN"_dna5)));
}
