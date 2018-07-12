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

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

#include <range/v3/view/take.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/range/view/translation.hpp>
#include <seqan3/std/concept/range.hpp>



using namespace seqan3;
using namespace seqan3::literal;
using namespace seqan3::view;

template <typename T>
class nucleotide : public ::testing::Test
{};

// add all alphabets here
using nucleotide_types  = ::testing::Types<dna4, dna5, dna15, rna4, rna5, rna15>;

TYPED_TEST_CASE(nucleotide, nucleotide_types);

TYPED_TEST(nucleotide, view_translate_single)
{
    using t = TypeParam;

    std::vector<TypeParam> vec = {t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A};
    aa27_vector cmp1{"TYVR"_aa27};
    aa27_vector cmp2{"CMHA"_aa27};

    // default parameter translation_frames
    auto v1 = vec | view::translate_single;
    // == [T,Y,V,R]
    EXPECT_TRUE((ranges::equal(v1, cmp1)));

    // default parameter translation_frames
//     auto v2 = vec | view::translate_single();
    // == [T,Y,V,R]
//     EXPECT_TRUE((ranges::equal(v2, cmp1));

    // single frame translation
    auto v3 = vec | view::translate_single(translation_frames::FWD_FRAME_0);
    // == [T,Y,V,R]
    EXPECT_TRUE((ranges::equal(v3, cmp1)));

    // function syntax
    auto v4 = view::translate_single(vec, translation_frames::FWD_FRAME_0);
    // == [T,Y,V,R]
    EXPECT_TRUE((ranges::equal(v4, cmp1)));

    // combinability
    auto v5 = vec | view::complement | view::translate_single(translation_frames::FWD_FRAME_0);
    // == [C,M,H,A]
    EXPECT_TRUE((ranges::equal(v5, cmp2)));
}

TYPED_TEST(nucleotide, view_translate)
{
    using t = TypeParam;

    std::vector<TypeParam> vec = {t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A};
    std::vector<std::vector<aa27> > cmp1{{"TYVR"_aa27}};
    std::vector<std::vector<aa27> > cmp2{{"TYVR"_aa27}, {"YVRT"_aa27}};
    std::vector<std::vector<aa27> > cmp3{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}};
    std::vector<std::vector<aa27> > cmp4{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}, {"YVRT"_aa27}, {"TYVR"_aa27}, {"RTY"_aa27}};
    std::vector<std::vector<aa27> > cmp5{{"TYVR"_aa27}, {"VRT"_aa27}};
    std::vector<std::vector<aa27> > cmp6{{"CMHA"_aa27}, {"MHAC"_aa27}};
    std::vector<std::vector<aa27> > cmp7{{"CMHA"_aa27}};

    // default parameter translation_frames
    auto v1 = vec | view::translate;
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]
    EXPECT_EQ(v1.size(), cmp4.size());
    for (unsigned i = 0; i < v1.size(); i++)
        EXPECT_TRUE((ranges::equal(v1[i], cmp4[i])));

    // default parameter translation_frames
//     auto v2 = vec | view::translate();
    // == [[T,Y,V,R]]
//     EXPECT_EQ(v2.size(), cmp4.size());
//     for (unsigned i = 0; i < v2.size(); i++)
//         EXPECT_TRUE((ranges::equal(v2[i], cmp4[i])));

    // single frame translation
    auto v3 = vec | view::translate(translation_frames::FWD_FRAME_0);
    // == [[T,Y,V,R]]
    EXPECT_EQ(v3.size(), cmp1.size());
    for (unsigned i = 0; i < v3.size(); i++)
        EXPECT_TRUE((ranges::equal(v3[i], cmp1[i])));

    // reverse translation
    auto v4 = vec | view::translate(translation_frames::FWD_REV_0);
    // == [[T,Y,V,R],[Y,V,R,T]]
    EXPECT_EQ(v4.size(), cmp2.size());
    for (unsigned i = 0; i < v4.size(); i++)
        EXPECT_TRUE((ranges::equal(v4[i], cmp2[i])));

    // forward frames translation
    auto v5 = vec | view::translate(translation_frames::FWD);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T]]
    EXPECT_EQ(v5.size(), cmp3.size());
    for (unsigned i = 0; i < v5.size(); i++)
        EXPECT_TRUE((ranges::equal(v5[i], cmp3[i])));

    // six frame translation
    auto v6 = vec | view::translate(translation_frames::SIX_FRAME);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]
    EXPECT_EQ(v6.size(), cmp4.size());
    for (unsigned i = 0; i < v6.size(); i++)
        EXPECT_TRUE((ranges::equal(v6[i], cmp4[i])));

    // user-defined frame combination
    auto v7 = vec | view::translate(translation_frames::FWD_FRAME_0 | translation_frames::FWD_FRAME_2);
    // == [[T,Y,V,R],[V,R,T]]
    EXPECT_EQ(v7.size(), cmp5.size());
    for (unsigned i = 0; i < v7.size(); i++)
        EXPECT_TRUE((ranges::equal(v7[i], cmp5[i])));

    // function syntax
    auto v8 = view::translate(vec, translation_frames::FWD_REV_0);
    // == [[T,Y,V,R],[Y,V,R,T]]
    EXPECT_EQ(v8.size(), cmp2.size());
    for (unsigned i = 0; i < v8.size(); i++)
        EXPECT_TRUE((ranges::equal(v8[i], cmp2[i])));

    // combinability
    auto v9 = vec | view::complement | view::translate(translation_frames::FWD_REV_0);
    // == [[C,M,H,A],[M,H,A,C]]
    EXPECT_EQ(v9.size(), cmp6.size());
    for (unsigned i = 0; i < v9.size(); i++)
        EXPECT_TRUE((ranges::equal(v9[i], cmp6[i])));

    // combinability
    auto v10 = vec | view::complement | view::translate(translation_frames::FWD_REV_0) | ranges::view::take(1);
    // == [[C,M,H,A]]
    EXPECT_EQ(v10.size(), cmp7.size());
    for (unsigned i = 0; i < v10.size(); i++)
        EXPECT_TRUE((ranges::equal(v10[i], cmp7[i])));

    // combinability and function syntax
    auto v11 = detail::view_translate(view::complement(vec), translation_frames::FWD_REV_0);
    // == [[T,Y,V,R],[Y,V,R,T]]
    EXPECT_EQ(v11.size(), cmp6.size());
    for (unsigned i = 0; i < v11.size(); i++)
        EXPECT_TRUE((ranges::equal(v11[i], cmp6[i])));
}

TYPED_TEST(nucleotide, view_translate_single_container_conversion)
{
    using t = TypeParam;

    std::vector<TypeParam> vec= {t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A};
    aa27_vector cmp1{"TYVR"_aa27};

    // default parameter translation_frames
    auto v1 = vec | view::translate_single;
    // == [T,Y,V,R]
    EXPECT_EQ(std::vector<aa27>(v1), cmp1);
}

TYPED_TEST(nucleotide, view_translate_container_conversion)
{
    using t = TypeParam;

    std::vector<TypeParam> vec= {t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A};
    std::vector<std::vector<aa27> > cmp1{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}, {"YVRT"_aa27}, {"TYVR"_aa27}, {"RTY"_aa27}};

    // six frame translation
    auto v1 = vec | view::translate(translation_frames::SIX_FRAME);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y]]
    EXPECT_EQ(v1.size(), cmp1.size());
    for (unsigned i = 0; i < v1.size(); i++)
        EXPECT_EQ(std::vector<std::vector<aa27> >(v1), cmp1);

    for (unsigned i = 0; i < v1.size(); i++)
        EXPECT_EQ(concatenated_sequences<std::vector<aa27> >(v1), cmp1);
}

TYPED_TEST(nucleotide, view_translate_single_concepts)
{
    using t = TypeParam;

    std::vector<TypeParam> vec= {t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A};
    EXPECT_TRUE(input_range_concept<decltype(vec)>);
    EXPECT_TRUE(forward_range_concept<decltype(vec)>);
    EXPECT_TRUE(random_access_range_concept<decltype(vec)>);
    EXPECT_TRUE(sized_range_concept<decltype(vec)>);

    auto v1 = vec | view::translate_single(translation_frames::FWD_FRAME_0);
    EXPECT_TRUE(input_range_concept<decltype(v1)>);
    EXPECT_TRUE(forward_range_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<decltype(v1)>);
    EXPECT_TRUE(sized_range_concept<decltype(v1)>);
    EXPECT_TRUE(view_concept<decltype(v1)>);
    EXPECT_TRUE((std::is_same_v<aa27, value_type_t<decltype(v1)>>));
    EXPECT_TRUE((std::is_same_v<aa27, reference_t<decltype(v1)>>));
}

TYPED_TEST(nucleotide, view_translate_concepts)
{
    using t = TypeParam;

    std::vector<TypeParam> vec= {t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A, t::C, t::G, t::T, t::A};
    EXPECT_TRUE(forward_range_concept<decltype(vec)>);
    EXPECT_TRUE(random_access_range_concept<decltype(vec)>);
    EXPECT_TRUE(sized_range_concept<decltype(vec)>);

    auto v1 = vec | view::translate(translation_frames::FWD_REV_0);
    EXPECT_TRUE(input_range_concept<decltype(v1)>);
    EXPECT_TRUE(forward_range_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<decltype(v1)>);
    EXPECT_TRUE(sized_range_concept<decltype(v1)>);
    EXPECT_TRUE(view_concept<decltype(v1)>);
    EXPECT_TRUE(random_access_range_concept<value_type_t<decltype(v1)>>);
    EXPECT_TRUE(sized_range_concept<value_type_t<decltype(v1)>>);
    EXPECT_TRUE(view_concept<value_type_t<decltype(v1)>>);
    EXPECT_TRUE(random_access_range_concept<reference_t<decltype(v1)>>);
    EXPECT_TRUE(sized_range_concept<reference_t<decltype(v1)>>);
    EXPECT_TRUE(view_concept<reference_t<decltype(v1)>>);
}
