// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate_join.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../../range/iterator_test_template.hpp"

using seqan3::operator""_aa27;
using seqan3::operator""_dna4;

using iterator_type =
    decltype(seqan3::views::translate_join(std::declval<std::vector<seqan3::dna4_vector> &>()).begin());

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<seqan3::dna4_vector> vec{"ACGTACGTACGTA"_dna4, "TCGAGAGCTTTAGC"_dna4};
    std::vector<std::vector<seqan3::aa27>> expected_range{{"TYVR"_aa27},
                                                          {"RTYV"_aa27},
                                                          {"VRT"_aa27},
                                                          {"YVRT"_aa27},
                                                          {"TYVR"_aa27},
                                                          {"RTY"_aa27},
                                                          {"SRAL"_aa27},
                                                          {"REL*"_aa27},
                                                          {"ESFS"_aa27},
                                                          {"AKAL"_aa27},
                                                          {"LKLS"_aa27},
                                                          {"*SSR"_aa27}};
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
using nucleotide_types =
    ::testing::Types<seqan3::dna4, seqan3::dna5, seqan3::dna15, seqan3::rna4, seqan3::rna5, seqan3::rna15>;

TYPED_TEST_SUITE(nucleotide, nucleotide_types, );

TYPED_TEST(nucleotide, view_translate)
{
    std::vector<std::string> sequences{"ACGTACGTACGTA", "TCGAGAGCTTTAGC"};
    auto vec = sequences | seqan3::views::char_to<TypeParam>;

    // default parameter translation_frames
    auto v1 = vec | seqan3::views::translate_join;
    EXPECT_EQ(v1.size(), 12u);
    EXPECT_RANGE_EQ("TYVR"_aa27, v1[0]);
    EXPECT_RANGE_EQ("RTYV"_aa27, v1[1]);
    EXPECT_RANGE_EQ("VRT"_aa27, v1[2]);
    EXPECT_RANGE_EQ("YVRT"_aa27, v1[3]);
    EXPECT_RANGE_EQ("TYVR"_aa27, v1[4]);
    EXPECT_RANGE_EQ("RTY"_aa27, v1[5]);
    EXPECT_RANGE_EQ("SRAL"_aa27, v1[6]);
    EXPECT_RANGE_EQ("REL*"_aa27, v1[7]);
    EXPECT_RANGE_EQ("ESFS"_aa27, v1[8]);
    EXPECT_RANGE_EQ("AKAL"_aa27, v1[9]);
    EXPECT_RANGE_EQ("LKLS"_aa27, v1[10]);
    EXPECT_RANGE_EQ("*SSR"_aa27, v1[11]);

    // default parameter translation_frames
    auto v2 = vec | seqan3::views::translate_join();
    EXPECT_EQ(v2.size(), 12u);
    EXPECT_RANGE_EQ("TYVR"_aa27, v2[0]);
    EXPECT_RANGE_EQ("RTYV"_aa27, v2[1]);
    EXPECT_RANGE_EQ("VRT"_aa27, v2[2]);
    EXPECT_RANGE_EQ("YVRT"_aa27, v2[3]);
    EXPECT_RANGE_EQ("TYVR"_aa27, v2[4]);
    EXPECT_RANGE_EQ("RTY"_aa27, v2[5]);
    EXPECT_RANGE_EQ("SRAL"_aa27, v2[6]);
    EXPECT_RANGE_EQ("REL*"_aa27, v2[7]);
    EXPECT_RANGE_EQ("ESFS"_aa27, v2[8]);
    EXPECT_RANGE_EQ("AKAL"_aa27, v2[9]);
    EXPECT_RANGE_EQ("LKLS"_aa27, v2[10]);
    EXPECT_RANGE_EQ("*SSR"_aa27, v2[11]);

    // single frame translation
    auto v3 = vec | seqan3::views::translate_join(seqan3::translation_frames::forward_frame0);
    EXPECT_EQ(v3.size(), 2u);
    EXPECT_RANGE_EQ(v3[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v3[1], "SRAL"_aa27);

    // reverse translation
    auto v4 = vec | seqan3::views::translate_join(seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v4.size(), 4u);
    EXPECT_RANGE_EQ(v4[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v4[1], "YVRT"_aa27);
    EXPECT_RANGE_EQ(v4[2], "SRAL"_aa27);
    EXPECT_RANGE_EQ(v4[3], "AKAL"_aa27);

    // forward frames translation
    auto v5 = vec | seqan3::views::translate_join(seqan3::translation_frames::forward_frames);
    EXPECT_EQ(v5.size(), 6u);
    EXPECT_RANGE_EQ(v5[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v5[1], "RTYV"_aa27);
    EXPECT_RANGE_EQ(v5[2], "VRT"_aa27);
    EXPECT_RANGE_EQ(v5[3], "SRAL"_aa27);
    EXPECT_RANGE_EQ(v5[4], "REL*"_aa27);
    EXPECT_RANGE_EQ(v5[5], "ESFS"_aa27);

    // six frame translation
    auto v6 = vec | seqan3::views::translate_join(seqan3::translation_frames::six_frames);
    EXPECT_EQ(v6.size(), 12u);
    EXPECT_RANGE_EQ("TYVR"_aa27, v6[0]);
    EXPECT_RANGE_EQ("RTYV"_aa27, v6[1]);
    EXPECT_RANGE_EQ("VRT"_aa27, v6[2]);
    EXPECT_RANGE_EQ("YVRT"_aa27, v6[3]);
    EXPECT_RANGE_EQ("TYVR"_aa27, v6[4]);
    EXPECT_RANGE_EQ("RTY"_aa27, v6[5]);
    EXPECT_RANGE_EQ("SRAL"_aa27, v6[6]);
    EXPECT_RANGE_EQ("REL*"_aa27, v6[7]);
    EXPECT_RANGE_EQ("ESFS"_aa27, v6[8]);
    EXPECT_RANGE_EQ("AKAL"_aa27, v6[9]);
    EXPECT_RANGE_EQ("LKLS"_aa27, v6[10]);
    EXPECT_RANGE_EQ("*SSR"_aa27, v6[11]);

    // user-defined frame combination
    auto v7 = vec
            | seqan3::views::translate_join(seqan3::translation_frames::forward_frame0
                                            | seqan3::translation_frames::forward_frame2);
    EXPECT_EQ(v7.size(), 4u);
    EXPECT_RANGE_EQ(v7[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v7[1], "VRT"_aa27);
    EXPECT_RANGE_EQ(v7[2], "SRAL"_aa27);
    EXPECT_RANGE_EQ(v7[3], "ESFS"_aa27);

    // function syntax
    auto v8 = seqan3::views::translate_join(vec, seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v8.size(), 4u);
    EXPECT_RANGE_EQ(v8[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v8[1], "YVRT"_aa27);
    EXPECT_RANGE_EQ(v8[2], "SRAL"_aa27);
    EXPECT_RANGE_EQ(v8[3], "AKAL"_aa27);

    // combinability
    auto v9 =
        vec | seqan3::views::complement | seqan3::views::translate_join(seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v9.size(), 4u);
    EXPECT_RANGE_EQ(v9[0], "CMHA"_aa27);
    EXPECT_RANGE_EQ(v9[1], "MHAC"_aa27);
    EXPECT_RANGE_EQ(v9[2], "SSRN"_aa27);
    EXPECT_RANGE_EQ(v9[3], "RFRE"_aa27);

    // combinability
    auto v10 = vec | seqan3::views::complement
             | seqan3::views::translate_join(seqan3::translation_frames::forward_reverse0) | std::views::take(1);
    EXPECT_EQ(v10.size(), 1u);
    EXPECT_RANGE_EQ(v10[0], "CMHA"_aa27);

    // combinability and function syntax
    auto v11 = seqan3::detail::view_translate_join(seqan3::views::complement(vec),
                                                   seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v11.size(), 4u);
    EXPECT_RANGE_EQ(v11[0], "CMHA"_aa27);
    EXPECT_RANGE_EQ(v11[1], "MHAC"_aa27);
    EXPECT_RANGE_EQ(v11[2], "SSRN"_aa27);
    EXPECT_RANGE_EQ(v11[3], "RFRE"_aa27);

    // combinability
    auto v12 = vec | seqan3::views::complement
             | seqan3::views::translate_join(seqan3::translation_frames::forward_reverse0)
             | seqan3::views::deep{std::views::reverse};
    EXPECT_EQ(v12.size(), 4u);
    EXPECT_RANGE_EQ(v12[0], "AHMC"_aa27);
    EXPECT_RANGE_EQ(v12[1], "CAHM"_aa27);
    EXPECT_RANGE_EQ(v12[2], "NRSS"_aa27);
    EXPECT_RANGE_EQ(v12[3], "ERFR"_aa27);
}

TYPED_TEST(nucleotide, view_translate_concepts)
{
    std::vector<std::vector<TypeParam>> vec{};

    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);

    auto v1 = vec | seqan3::views::translate_join(seqan3::translation_frames::forward_reverse0);

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
