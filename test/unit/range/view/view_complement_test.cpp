// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;
using namespace seqan3::literal;

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
    dna5_vector v3 = foo | view::complement | view::reverse;
    EXPECT_EQ(v3, "TACGT"_dna5);
}

TEST(view_complement, deep_view)
{
    std::vector<dna5_vector> foo{"ACGTA"_dna5, "TGCAT"_dna5};

    auto v = foo | view::complement;

    ASSERT_EQ(ranges::size(v), 2);
    EXPECT_TRUE((ranges::equal(v[0], "TGCAT"_dna5)));
    EXPECT_TRUE((ranges::equal(v[1], "ACGTA"_dna5)));
}

TEST(view_complement, concepts)
{
    dna5_vector vec{"ACGTA"_dna5};
    EXPECT_TRUE(input_range_concept<decltype(vec)>);
    EXPECT_TRUE(forward_range_concept<decltype(vec)>);
    EXPECT_TRUE(bidirectional_range_concept<decltype(vec)>);
    EXPECT_TRUE(random_access_range_concept<decltype(vec)>);
    EXPECT_FALSE(view_concept<decltype(vec)>);
    EXPECT_TRUE(sized_range_concept<decltype(vec)>);
    EXPECT_TRUE(common_range_concept<decltype(vec)>);
    EXPECT_TRUE(const_iterable_concept<decltype(vec)>);
    EXPECT_TRUE((output_range_concept<decltype(vec), dna5>));

    auto v1 = vec | view::complement;
    EXPECT_TRUE(input_range_concept<decltype(v1)>);
    EXPECT_TRUE(forward_range_concept<decltype(v1)>);
    EXPECT_TRUE(bidirectional_range_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<decltype(v1)>);
    EXPECT_TRUE(view_concept<decltype(v1)>);
    EXPECT_TRUE(sized_range_concept<decltype(v1)>);
    EXPECT_TRUE(common_range_concept<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_FALSE((output_range_concept<decltype(v1), dna5>));
    EXPECT_FALSE((output_range_concept<decltype(v1), char>));
}
