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

#include <range/v3/view/reverse.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/trim.hpp>

using namespace seqan3;
using namespace seqan3::view;

TEST(view_trim, standalone)
{
    std::vector<illumina18> vec{ illumina18{40}, illumina18{40}, illumina18{30}, illumina18{20}, illumina18{10}};
    std::vector<illumina18> cmp1{illumina18{40}, illumina18{40}, illumina18{30}, illumina18{20}};
    std::vector<illumina18> cmp2{illumina18{40}, illumina18{40}};

    // trim by phred_value
    auto v1 = vec | view::trim(20u);                        // == ['I','I','?','5']
    EXPECT_EQ(std::vector<illumina18>(v1), cmp1);

    // trim by quality character
    auto v2 = vec | view::trim(illumina18{40});             // == ['I','I']
    EXPECT_EQ(std::vector<illumina18>(v2), cmp2);

    // function syntax
    auto v3 = view::trim(vec, 20u);                         // == ['I','I','?','5']
    EXPECT_EQ(std::vector<illumina18>(v3), cmp1);

    // combinability
    std::string v4 = view::trim(vec, 20u) | view::to_char;  // == "II?5"
    EXPECT_EQ("II?5", v4);
}

TEST(view_trim, quality_composition)
{
    std::vector<dna5q> vec{{dna5::A, illumina18{40}}, {dna5::G, illumina18{40}}, {dna5::G, illumina18{30}},
                           {dna5::A, illumina18{20}}, {dna5::T, illumina18{10}}};
    std::vector<dna5q> cmp1{{dna5::A, illumina18{40}}, {dna5::G, illumina18{40}}, {dna5::G, illumina18{30}},
                            {dna5::A, illumina18{20}}};
    std::vector<dna5q> cmp2{{dna5::A, illumina18{40}}, {dna5::G, illumina18{40}}};

    // trim by phred_value
    auto v1 = vec | view::trim(20u);
    EXPECT_EQ(std::vector<dna5q>(v1), cmp1);

    // trim by quality character
    auto v2 = vec | view::trim(dna5q{dna5::C, 40});
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
    std::vector<dna5q> vec{{dna5::A, 40}, {dna5::G, 40}, {dna5::G, 30}, {dna5::A, 20}, {dna5::T, 10}};
    EXPECT_TRUE(input_range_concept<decltype(vec)>);
    EXPECT_TRUE(forward_range_concept<decltype(vec)>);
    EXPECT_TRUE(random_access_range_concept<decltype(vec)>);
    EXPECT_TRUE(sized_range_concept<decltype(vec)>);

    auto v1 = vec | view::trim(20u);
    EXPECT_TRUE(input_range_concept<decltype(v1)>);
    EXPECT_TRUE(forward_range_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<decltype(v1)>);
    EXPECT_TRUE(!sized_range_concept<decltype(v1)>);
}
