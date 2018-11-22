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

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/trim.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;
using namespace seqan3::view;

TEST(view_trim, standalone)
{
    std::vector<phred42> vec{ phred42{40}, phred42{40}, phred42{30}, phred42{20}, phred42{10}};
    std::vector<phred42> cmp1{phred42{40}, phred42{40}, phred42{30}, phred42{20}};
    std::vector<phred42> cmp2{phred42{40}, phred42{40}};

    // trim by phred_value
    auto v1 = vec | view::trim(20u);                        // == ['I','I','?','5']
    EXPECT_EQ(std::vector<phred42>(v1), cmp1);

    // trim by quality character
    auto v2 = vec | view::trim(phred42{40});             // == ['I','I']
    EXPECT_EQ(std::vector<phred42>(v2), cmp2);

    // function syntax
    auto v3 = view::trim(vec, 20u);                         // == ['I','I','?','5']
    EXPECT_EQ(std::vector<phred42>(v3), cmp1);

    // combinability
    std::string v4 = view::trim(vec, 20u) | view::to_char;  // == "II?5"
    EXPECT_EQ("II?5", v4);
}

TEST(view_trim, qualified)
{
    std::vector<dna5q> vec{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}},
                           {'A'_dna5, phred42{20}}, {'T'_dna5, phred42{10}}};
    std::vector<dna5q> cmp1{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}},
                            {'A'_dna5, phred42{20}}};
    std::vector<dna5q> cmp2{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}};

    // trim by phred_value
    auto v1 = vec | view::trim(20u);
    EXPECT_EQ(std::vector<dna5q>(v1), cmp1);

    // trim by quality character
    auto v2 = vec | view::trim(dna5q{'C'_dna5, phred42{40}});
    EXPECT_EQ(std::vector<dna5q>(v2), cmp2);

    // function syntax
    auto v3 = view::trim(vec, 20u);
    EXPECT_EQ(std::vector<dna5q>(v3), cmp1);

    // combinability
    std::string v4 = view::trim(vec, 20u) | view::to_char;
    EXPECT_EQ("AGGA", v4);
}

TEST(view_trim, concepts)
{
    std::vector<dna5q> vec{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}}, {'A'_dna5, phred42{20}}, {'T'_dna5, phred42{10}}};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), dna5q>));
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);

    auto v1 = vec | view::trim(20u);
    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), dna5q>));
    EXPECT_TRUE(!std::ranges::SizedRange<decltype(v1)>);
}
