// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>

#include <gtest/gtest.h>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

static constexpr auto ungapped_view = seqan3::views::kmer_hash(seqan3::ungapped{3});
static constexpr auto gapped_view = seqan3::views::kmer_hash(0b101_shape);
static constexpr auto prefix_until_first_thymine = seqan3::views::take_until([] (seqan3::dna4 x)
                                                   { return x == 'T'_dna4; });

template <typename T>
class kmer_hash_test: public ::testing::Test {};

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitcompressed_vector<seqan3::dna4>,
                                                seqan3::bitcompressed_vector<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;

TYPED_TEST_SUITE(kmer_hash_test, underlying_range_types, );

TYPED_TEST(kmer_hash_test, ungapped_combined_with_container)
{
    TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
    result_t ungapped1{6, 27, 44, 50, 9};
    TypeParam text2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4}; // AAAAA
    result_t ungapped2{0,0,0};
    TypeParam text3{'A'_dna4, 'C'_dna4}; // AC
    result_t ungapped3{};
    TypeParam text4{'A'_dna4, 'C'_dna4, 'G'_dna4}; // AC
    result_t ungapped4{6};

    EXPECT_EQ(ungapped1, text1 | ungapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped2, text2 | ungapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped3, text3 | ungapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped4, text4 | ungapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped4, text1 | prefix_until_first_thymine | ungapped_view | seqan3::views::to<result_t>);
}

TYPED_TEST(kmer_hash_test, gapped_combined_with_container)
{
    TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
    result_t gapped1{2, 7, 8, 14, 1};
    TypeParam text2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4}; // AAAAA
    result_t gapped2{0,0,0};
    TypeParam text3{'A'_dna4, 'C'_dna4}; // AC
    result_t gapped3{};
    TypeParam text4{'A'_dna4, 'C'_dna4, 'G'_dna4}; // AC
    result_t gapped4{2};
    EXPECT_EQ(gapped1, text1 | gapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped2, text2 | gapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped3, text3 | gapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped4, text4 | gapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped4, text1 | prefix_until_first_thymine| gapped_view | seqan3::views::to<result_t>);
}

TYPED_TEST(kmer_hash_test, ungapped_concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}; // ACGT
    auto v1 = text | ungapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::bidirectional_range<decltype(text)>, std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(text)>, std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));
}

TYPED_TEST(kmer_hash_test, gapped_concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}; // ACGT
    auto v1 = text | gapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::bidirectional_range<decltype(text)>, std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(text)>, std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));
}

TYPED_TEST(kmer_hash_test, invalid_sizes)
{
    TypeParam text1{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    EXPECT_NO_THROW(text1 | seqan3::views::kmer_hash(seqan3::ungapped{32}));
    EXPECT_THROW(text1 | seqan3::views::kmer_hash(seqan3::ungapped{33}), std::invalid_argument);
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
       EXPECT_NO_THROW(text1 | std::views::reverse | seqan3::views::kmer_hash(seqan3::ungapped{32}));
       EXPECT_THROW(text1 | std::views::reverse | seqan3::views::kmer_hash(seqan3::ungapped{33}),
                    std::invalid_argument);
    }

    EXPECT_NO_THROW(text1 | seqan3::views::kmer_hash(0xFFFFFFFE001_shape)); // size=44, count=32
    EXPECT_THROW(text1 | seqan3::views::kmer_hash(0xFFFFFFFFE009_shape), std::invalid_argument); // size=44, count=33

    std::vector<seqan3::dna5> dna5_text{};
    EXPECT_NO_THROW(dna5_text | seqan3::views::kmer_hash(seqan3::ungapped{27}));
    EXPECT_THROW(dna5_text | seqan3::views::kmer_hash(seqan3::ungapped{28}), std::invalid_argument);
}

// https://github.com/seqan/seqan3/issues/1614
TEST(kmer_hash_test, issue1614)
{
    std::vector<seqan3::dna5> sequence{"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5};
    EXPECT_EQ(sequence | seqan3::views::kmer_hash(seqan3::ungapped{25}) | seqan3::views::to<std::vector<size_t>>,
              seqan3::views::repeat_n(298023223876953124, 26) | seqan3::views::to<std::vector<size_t>>);
}

// https://github.com/seqan/seqan3/issues/1643
TEST(kmer_hash_test, issue1643)
{
    std::vector<seqan3::dna4> text_23_elements{"ACGATCGATCGTAGCTACTGAGC"_dna4};

    auto k_mer_size_23_view = text_23_elements | seqan3::views::kmer_hash(seqan3::ungapped{23u});
    EXPECT_EQ(k_mer_size_23_view.size(), 1u);
    EXPECT_EQ(k_mer_size_23_view[0], 6829917194121u);

    auto k_mer_size_24_view = text_23_elements | seqan3::views::kmer_hash(seqan3::ungapped{24u});
    EXPECT_TRUE(k_mer_size_24_view.empty());

    auto k_mer_size_25_view = text_23_elements | seqan3::views::kmer_hash(seqan3::ungapped{25u});
    EXPECT_TRUE(k_mer_size_25_view.empty());
}

// https://github.com/seqan/seqan3/issues/1719
TEST(kmer_hash_test, issue1719)
{
    uint64_t const expected = 0;
    std::vector<seqan3::dna5> sequence{""_dna5};
    auto v = sequence | seqan3::views::kmer_hash(seqan3::ungapped{25});
    EXPECT_EQ(expected, v.size());

    std::vector<seqan3::dna5> sequence2{"ACGATCGATCGTAGCTACTGAGC"_dna5};
    auto v2 = sequence2 | seqan3::views::kmer_hash(seqan3::ungapped{25});
    EXPECT_EQ(expected, v2.size());
}
