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
#include <seqan3/test/expect_range_eq.hpp>

#include "../iterator_test_template.hpp"

#include <gtest/gtest.h>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

static constexpr auto ungapped_view = seqan3::views::kmer_hash(seqan3::ungapped{3});
static constexpr auto gapped_view = seqan3::views::kmer_hash(0b101_shape);
static constexpr auto prefix_until_first_thymine = seqan3::views::take_until([] (seqan3::dna4 x)
                                                   { return x == 'T'_dna4; });

using iterator_type = std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector&>() | gapped_view)>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGTAGC"_dna4};

    decltype(text | gapped_view) test_range = text | gapped_view;

    std::vector<size_t> expected_range{2, 7, 8, 14, 1};
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class kmer_hash_ungapped_test: public ::testing::Test {};

template <typename T>
class kmer_hash_gapped_test: public ::testing::Test {};

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitcompressed_vector<seqan3::dna4>,
                                                seqan3::bitcompressed_vector<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;

TYPED_TEST_SUITE(kmer_hash_ungapped_test, underlying_range_types, );
TYPED_TEST_SUITE(kmer_hash_gapped_test, underlying_range_types, );

TYPED_TEST(kmer_hash_ungapped_test, combined_with_container)
{
    {
        TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
        result_t ungapped1{6, 27, 44, 50, 9};
        EXPECT_RANGE_EQ(ungapped1, text1 | ungapped_view);
        EXPECT_RANGE_EQ(result_t{6}, text1 | prefix_until_first_thymine | ungapped_view);
    }
    {
        TypeParam text2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4}; // AAAAA
        result_t ungapped2{0, 0, 0};
        EXPECT_RANGE_EQ(ungapped2, text2 | ungapped_view);
    }
    {
        TypeParam text3{'A'_dna4, 'C'_dna4}; // AC
        EXPECT_RANGE_EQ(result_t{}, text3 | ungapped_view);
    }
    {
        TypeParam text4{'A'_dna4, 'C'_dna4, 'G'_dna4}; // AC
        EXPECT_RANGE_EQ(result_t{6}, text4 | ungapped_view);
    }
}

TYPED_TEST(kmer_hash_gapped_test, combined_with_container)
{
    {
        TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
        result_t gapped1{2, 7, 8, 14, 1};
        EXPECT_RANGE_EQ(gapped1, text1 | gapped_view);
        EXPECT_RANGE_EQ(result_t{2}, text1 | prefix_until_first_thymine| gapped_view);
    }
    {
        TypeParam text2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4}; // AAAAA
        result_t gapped2{0, 0, 0};
        EXPECT_RANGE_EQ(gapped2, text2 | gapped_view);
    }
    {
        TypeParam text3{'A'_dna4, 'C'_dna4}; // AC
        EXPECT_RANGE_EQ(result_t{}, text3 | gapped_view);
    }
    {
        TypeParam text4{'A'_dna4, 'C'_dna4, 'G'_dna4}; // AC
        EXPECT_RANGE_EQ(result_t{2}, text4 | gapped_view);
    }
}

TYPED_TEST(kmer_hash_ungapped_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}; // ACGT
    auto v1 = text | ungapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::bidirectional_range<decltype(text)>, std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(text)>, std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));
}

TYPED_TEST(kmer_hash_gapped_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4}; // ACGT
    auto v1 = text | gapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::bidirectional_range<decltype(text)>, std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(text)>, std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(text)>, std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));
}

TYPED_TEST(kmer_hash_ungapped_test, invalid_sizes)
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
TEST(kmer_hash_ungapped_test, issue1614)
{
    std::vector<seqan3::dna5> sequence{"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5};
    EXPECT_RANGE_EQ(sequence | seqan3::views::kmer_hash(seqan3::ungapped{25}),
                    seqan3::views::repeat_n(298023223876953124, 26));
}

// https://github.com/seqan/seqan3/issues/1643
TEST(kmer_hash_ungapped_test, issue1643)
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
TYPED_TEST(kmer_hash_ungapped_test, issue1719)
{
    if constexpr (std::ranges::sized_range<TypeParam>)
    {
        TypeParam sequence{};
        auto v = sequence | seqan3::views::kmer_hash(seqan3::ungapped{8});
        EXPECT_EQ(0u, v.size());

        TypeParam sequence2{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4};
        auto v2 = sequence2 | seqan3::views::kmer_hash(seqan3::ungapped{8});
        EXPECT_EQ(0u, v2.size());

        auto v3 = sequence2 | seqan3::views::kmer_hash(seqan3::ungapped{4});
        EXPECT_EQ(4u, v3.size());
    }
}
