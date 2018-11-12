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

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;
using namespace seqan3::literal;

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
    std::vector<bool> v3 = vec | view::convert<bool> | view::reverse;
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
    dna4_vector v3 = vec | view::convert<dna4> | view::reverse;
    EXPECT_EQ(cmp2, v3);
}

TEST(view_convert, concepts)
{
    dna5_vector vec{"ACGNTNGGN"_dna5};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(ConstIterableRange<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), dna5>));

    auto v1 = vec | view::convert<dna4>;
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(ConstIterableRange<decltype(v1)>);
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna5>));
    EXPECT_FALSE((std::ranges::OutputRange<decltype(v1), dna4>));
}
