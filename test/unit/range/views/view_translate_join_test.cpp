// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/translate_join.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../iterator_test_template.hpp"

using seqan3::operator""_aa27;
using seqan3::operator""_dna4;

using iterator_type =
        decltype(seqan3::views::translate_join(std::declval<std::vector<seqan3::dna4_vector > &>()).begin());

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<seqan3::dna4_vector > vec{"ACGTACGTACGTA"_dna4, "TCGAGAGCTTTAGC"_dna4};
    std::vector<std::vector<seqan3::aa27> > expected_range{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}, {"YVRT"_aa27},
                                                           {"TYVR"_aa27}, {"RTY"_aa27}, {"SRAL"_aa27}, {"REL*"_aa27},
                                                           {"ESFS"_aa27}, {"AKAL"_aa27}, {"LKLS"_aa27}, {"*SSR"_aa27}};
    decltype(seqan3::views::translate_join(vec)) test_range = seqan3::views::translate_join(vec);

    template <typename A, typename B>
    static void expect_eq(A && test_range_value, B && expected_range_value)
    {
        EXPECT_RANGE_EQ(test_range_value, expected_range_value);
    }
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class nucleotide : public ::testing::Test
{};

// add all alphabets here
using nucleotide_types = ::testing::Types<seqan3::dna4,
                                          seqan3::dna5,
                                          seqan3::dna15,
                                          seqan3::rna4,
                                          seqan3::rna5,
                                          seqan3::rna15>;

TYPED_TEST_SUITE(nucleotide, nucleotide_types, );

TYPED_TEST(nucleotide, view_translate)
{
    std::string const in1{"ACGTACGTACGTA"};
    std::string const in2{"TCGAGAGCTTTAGC"};
    std::vector<std::vector<TypeParam> > vec;
    vec.resize(2);
    vec[0] = in1 | seqan3::views::char_to<TypeParam> | seqan3::views::to<std::vector>;
    vec[1] = in2 | seqan3::views::char_to<TypeParam> | seqan3::views::to<std::vector>;

    std::vector<std::vector<seqan3::aa27> > cmp1{{"TYVR"_aa27}, {"SRAL"_aa27}};
    std::vector<std::vector<seqan3::aa27> > cmp2{{"TYVR"_aa27}, {"YVRT"_aa27}, {"SRAL"_aa27}, {"AKAL"_aa27}};
    std::vector<std::vector<seqan3::aa27> > cmp3{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}, {"SRAL"_aa27},
                                                 {"REL*"_aa27}, {"ESFS"_aa27}};
    std::vector<std::vector<seqan3::aa27> > cmp4{{"TYVR"_aa27}, {"RTYV"_aa27}, {"VRT"_aa27}, {"YVRT"_aa27},
                                                 {"TYVR"_aa27}, {"RTY"_aa27}, {"SRAL"_aa27}, {"REL*"_aa27},
                                                 {"ESFS"_aa27}, {"AKAL"_aa27}, {"LKLS"_aa27}, {"*SSR"_aa27}};
    std::vector<std::vector<seqan3::aa27> > cmp5{{"TYVR"_aa27}, {"VRT"_aa27}, {"SRAL"_aa27}, {"ESFS"_aa27}};
    std::vector<std::vector<seqan3::aa27> > cmp6{{"CMHA"_aa27}, {"MHAC"_aa27}, {"SSRN"_aa27}, {"RFRE"_aa27}};
    std::vector<std::vector<seqan3::aa27> > cmp7{{"CMHA"_aa27}};

    // default parameter translation_frames
    auto v1 = vec | seqan3::views::translate_join;
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y],[S,R,A,L],[R,E,L,*],[E,S,F,S],[A,K,A,L],[L,K,L,S],[*,S,S,R]]
    EXPECT_EQ(v1.size(), cmp4.size());
    for (unsigned i = 0; i < v1.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v1[i], cmp4[i])));

    // default parameter translation_frames
    auto v2 = vec | seqan3::views::translate_join();
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y],[S,R,A,L],[R,E,L,*],[E,S,F,S],[A,K,A,L],[L,K,L,S],[*,S,S,R]]
    EXPECT_EQ(v2.size(), cmp4.size());
    for (unsigned i = 0; i < v2.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v2[i], cmp4[i])));

    // single frame translation
    auto v3 = vec | seqan3::views::translate_join(seqan3::translation_frames::FWD_FRAME_0);
    // == [[T,Y,V,R],[S,R,A,L]]
    EXPECT_EQ(v3.size(), cmp1.size());
    for (unsigned i = 0; i < v3.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v3[i], cmp1[i])));

    // reverse translation
    auto v4 = vec | seqan3::views::translate_join(seqan3::translation_frames::FWD_REV_0);
    // == [[T,Y,V,R],Y,V,R,T],[S,R,A,L],[A,K,A,L]]
    EXPECT_EQ(v4.size(), cmp2.size());
    for (unsigned i = 0; i < v4.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v4[i], cmp2[i])));

    // forward frames translation
    auto v5 = vec | seqan3::views::translate_join(seqan3::translation_frames::FWD);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[S,R,A,L],[R,E,L,*],[E,S,F,S]]
    EXPECT_EQ(v5.size(), cmp3.size());
    for (unsigned i = 0; i < v5.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v5[i], cmp3[i])));

    // six frame translation
    auto v6 = vec | seqan3::views::translate_join(seqan3::translation_frames::SIX_FRAME);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y],[S,R,A,L],[R,E,L,*],[E,S,F,S],[A,K,A,L],[L,K,L,S],[*,S,S,R]]
    EXPECT_EQ(v6.size(), cmp4.size());
    for (unsigned i = 0; i < v6.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v6[i], cmp4[i])));

    // user-defined frame combination
    auto v7 = vec
            | seqan3::views::translate_join(seqan3::translation_frames::FWD_FRAME_0
            | seqan3::translation_frames::FWD_FRAME_2);
    // == [[T,Y,V,R],[V,R,T],[S,R,A,L],[E,S,F,S]]
    EXPECT_EQ(v7.size(), cmp5.size());
    for (unsigned i = 0; i < v7.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v7[i], cmp5[i])));

    // function syntax
    auto v8 = seqan3::views::translate_join(vec, seqan3::translation_frames::FWD_REV_0);
    // == [[T,Y,V,R],[R,T,Y,V],[V,R,T],[Y,V,R,T],[T,Y,V,R],[R,T,Y],[S,R,A,L], [R,E,L,*], [E,S,F,S], [A,K,A,L], [L,K,L,S], [*,S,S,R]]
    EXPECT_EQ(v8.size(), cmp2.size());
    for (unsigned i = 0; i < v8.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v8[i], cmp2[i])));

    // combinability
    auto v9 = vec | seqan3::views::complement | seqan3::views::translate_join(seqan3::translation_frames::FWD_REV_0);
    // == [[C,M,H,A],[M,H,A,C],[S,S,R,N],[R,F,R,E]]
    EXPECT_EQ(v9.size(), cmp6.size());
    for (unsigned i = 0; i < v9.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v9[i], cmp6[i])));

    // combinability
    auto v10 = vec
             | seqan3::views::complement
             | seqan3::views::translate_join(seqan3::translation_frames::FWD_REV_0)
             | std::views::take(1);
    // == [C,M,H,A]
    EXPECT_EQ(v10.size(), cmp7.size());
    for (unsigned i = 0; i < v10.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v10[i], cmp7[i])));

    // combinability and function syntax
    auto v11 = seqan3::detail::view_translate_join(seqan3::views::complement(vec),
                                                   seqan3::translation_frames::FWD_REV_0);
    // == [[C,M,H,A],[M,H,A,C],[S,S,R,N],[R,F,R,E]]
    EXPECT_EQ(v11.size(), cmp6.size());
    for (unsigned i = 0; i < v11.size(); i++)
        EXPECT_TRUE((std::ranges::equal(v11[i], cmp6[i])));
}

TYPED_TEST(nucleotide, view_translate_concepts)
{
    std::string const in1{"ACGTACGTACGTA"};
    std::string const in2{"TCGAGAGCTTTAGC"};
    std::vector<std::vector<TypeParam> > vec;
    vec.resize(2);
    vec[0] = in1 | seqan3::views::char_to<TypeParam> | seqan3::views::to<std::vector>;
    vec[1] = in2 | seqan3::views::char_to<TypeParam> | seqan3::views::to<std::vector>;

    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);

    auto v1 = vec | seqan3::views::translate_join(seqan3::translation_frames::FWD_REV_0);

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<std::ranges::range_value_t<decltype(v1)>>);
    EXPECT_TRUE(std::ranges::sized_range<std::ranges::range_value_t<decltype(v1)>>);
    EXPECT_TRUE(std::ranges::view<std::ranges::range_value_t<decltype(v1)>>);
    EXPECT_TRUE(std::ranges::random_access_range<std::ranges::range_reference_t<decltype(v1)>>);
    EXPECT_TRUE(std::ranges::sized_range<std::ranges::range_reference_t<decltype(v1)>>);
    EXPECT_TRUE(std::ranges::view<std::ranges::range_reference_t<decltype(v1)>>);
}
