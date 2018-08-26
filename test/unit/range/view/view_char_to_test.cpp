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

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;
using namespace seqan3::literal;

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
    dna5_vector v3 = vec | view::char_to<dna5> | view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_char_to, deep_view)
{
    std::vector<std::string> foo{"ACGTA", "TGCAT"};

    std::vector<dna5_vector> v = foo | view::char_to<dna5>;

    ASSERT_EQ(ranges::size(v), 2);
    EXPECT_TRUE((ranges::equal(v[0], "ACGTA"_dna5)));
    EXPECT_TRUE((ranges::equal(v[1], "TGCAT"_dna5)));
}

TEST(view_char_to, concepts)
{
    std::string vec{"ACTTTGATA"};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), char>));

    auto v1 = vec | view::char_to<dna5>;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), char>));
}
