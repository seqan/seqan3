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
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/range/view/concept.hpp>
#include <seqan3/range/view/translation.hpp>

using namespace seqan3;
using namespace seqan3::literal;
using namespace seqan3::view;

TEST(view_translate_single, standalone)
{
    dna5_vector vec{"ACGTACGTACGTA"_dna5};
    aa27_vector cmp1{"TYVR"_aa27};
    aa27_vector cmp2{"CMHA"_aa27};

    // single frame translation
    auto v1 = vec | view::translate_single;                     // == [T,Y,V,R]
    EXPECT_EQ(std::vector<aa27>(v1), cmp1);

    // function syntax
    auto v2 = view::translate_single(vec);                      // == [T,Y,V,R]
    EXPECT_EQ(std::vector<aa27>(v2), cmp1);

    // combinability
    auto v3 = vec | view::complement | view::translate_single;  // == [C,M,H,A]
    EXPECT_EQ(std::vector<aa27>(v3), cmp2);
}

TEST(view_translate_frames, standalone)
{
    dna5_vector vec{"ACGTACGTACGTA"_dna5};
    std::vector<std::vector<aa27> > cmp1{{"TYVR"_aa27}};
    std::vector<std::vector<aa27> > cmp2{{"TYVR"_aa27}, {"YVRT"_aa27}};
    std::vector<std::vector<aa27> > cmp3{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}};
    std::vector<std::vector<aa27> > cmp4{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}, {"YVRT"_aa27}, {"TYVR"_aa27}, {"RTY"_aa27}};
    std::vector<std::vector<aa27> > cmp5{{"CMHA"_aa27}, {"MHAC"_aa27}};

    // single frame translation
    auto v1 = vec | view::translate_frames(translation_frames::SINGLE_FRAME);                                              // == [[T,Y,V,R]]
    EXPECT_EQ(std::vector<std::vector<aa27> >(v1), cmp1);

    // reverse translation
    auto v2 = vec | view::translate_frames(translation_frames::WITH_REVERSE_COMPLEMENT);                         // == [[T,Y,V,R],[Y,V,R,T]]
    EXPECT_EQ(std::vector<std::vector<aa27> >(v2), cmp2);

    // forward frames translation
    auto v3 = vec | view::translate_frames(translation_frames::WITH_FRAME_SHIFTS);                       // == [[T,Y,V,R],[R,T,Y,V],[V,R,T]]
    EXPECT_EQ(std::vector<std::vector<aa27> >(v3), cmp3);

    // six frame translation
    auto v4 = vec | view::translate_frames(translation_frames::SIX_FRAME);   // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]
    EXPECT_EQ(std::vector<std::vector<aa27> >(v4), cmp4);

    // function syntax
    auto v5 = view::translate_frames(vec, translation_frames::WITH_REVERSE_COMPLEMENT);                          // == [[T,Y,V,R],[Y,V,R,T]]
    EXPECT_EQ(std::vector<std::vector<aa27> >(v5), cmp2);

    // combinability
    auto v6 = vec | view::complement | view::translate_frames(translation_frames::WITH_REVERSE_COMPLEMENT);      // == [[C,M,H,A],[M,H,A,C]]
    EXPECT_EQ(std::vector<std::vector<aa27> >(v6), cmp5);
}

TEST(view_translate_single, concepts)
{
    dna5_vector vec{"ACGTACGTACGTA"_dna5};
    EXPECT_TRUE(input_range_concept<decltype(vec)>);
    EXPECT_TRUE(forward_range_concept<decltype(vec)>);
    EXPECT_TRUE(random_access_range_concept<decltype(vec)>);
    EXPECT_TRUE(sized_range_concept<decltype(vec)>);

    auto v1 = vec | view::translate_single;
    EXPECT_TRUE(input_range_concept<decltype(v1)>);
    EXPECT_TRUE(forward_range_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<decltype(v1)>);
    EXPECT_TRUE(sized_range_concept<decltype(v1)>);
    EXPECT_TRUE(view_concept<decltype(v1)>);
}

TEST(view_translate_frames, concepts)
{
    dna5_vector vec{"ACGTACGTACGTA"_dna5};
    EXPECT_TRUE(forward_range_concept<decltype(vec)>);
    EXPECT_TRUE(random_access_range_concept<decltype(vec)>);
    EXPECT_TRUE(sized_range_concept<decltype(vec)>);

    auto v1 = vec | view::translate_frames(translation_frames::WITH_REVERSE_COMPLEMENT);
    EXPECT_TRUE(input_range_concept<decltype(v1)>);
    EXPECT_TRUE(forward_range_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<decltype(v1)>);
    EXPECT_TRUE(!sized_range_concept<decltype(v1)>);
    EXPECT_TRUE(view_concept<decltype(v1)>);
}
