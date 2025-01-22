// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_throw_msg.hpp>
#include <seqan3/utility/views/repeat_n.hpp>

#include "../../range/iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

static constexpr auto ungapped_view = seqan3::views::kmer_hash(seqan3::ungapped{3});
static constexpr auto gapped_view = seqan3::views::kmer_hash(0b101_shape);
static constexpr auto prefix_until_first_thymine = std::views::take_while(
    [](seqan3::dna4 const x)
    {
        return x != 'T'_dna4;
    });

using iterator_type = std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector &>() | gapped_view)>;

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
class kmer_hash_ungapped_test : public ::testing::Test
{};

template <typename T>
class kmer_hash_gapped_test : public ::testing::Test
{};

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
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
        EXPECT_RANGE_EQ(result_t{2}, text1 | prefix_until_first_thymine | gapped_view);
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
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
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
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));
}

TYPED_TEST(kmer_hash_ungapped_test, invalid_sizes)
{
    auto expected_error_message =
        [](std::string_view const alphabet, size_t const max_shape_count, size_t const given_shape_count)
    {
        std::string message{"The shape is too long for the given alphabet.\n"};
        message += "Alphabet: ";
        message += alphabet;
        message += "\nMaximum shape count: ";
        message += std::to_string(max_shape_count);
        message += "\nGiven shape count: ";
        message += std::to_string(given_shape_count);
        return message;
    };

    TypeParam text{};
    EXPECT_NO_THROW(text | seqan3::views::kmer_hash(seqan3::ungapped{32}));
    EXPECT_THROW_MSG(text | seqan3::views::kmer_hash(seqan3::ungapped{33}),
                     std::invalid_argument,
                     expected_error_message("seqan3::dna4", 32, 33));

    if constexpr (std::ranges::bidirectional_range<TypeParam>)
    {
        EXPECT_NO_THROW(text | std::views::reverse | seqan3::views::kmer_hash(seqan3::ungapped{32}));
        EXPECT_THROW_MSG(text | std::views::reverse | seqan3::views::kmer_hash(seqan3::ungapped{33}),
                         std::invalid_argument,
                         expected_error_message("seqan3::dna4", 32, 33));
    }

    EXPECT_NO_THROW(text | seqan3::views::kmer_hash(0xF'FF'FF'FF'E0'01_shape));  // size=48, count=32
    EXPECT_THROW_MSG(text | seqan3::views::kmer_hash(0xFF'FF'FF'FE'00'09_shape), // size=48, count=33
                     std::invalid_argument,
                     expected_error_message("seqan3::dna4", 32, 33));

    std::vector<seqan3::dna5> dna5_text{};
    EXPECT_NO_THROW(dna5_text | seqan3::views::kmer_hash(seqan3::ungapped{27}));
    EXPECT_THROW_MSG(dna5_text | seqan3::views::kmer_hash(seqan3::ungapped{28}),
                     std::invalid_argument,
                     expected_error_message("seqan3::dna5", 27, 28));
}

// https://github.com/seqan/seqan3/issues/1614
TEST(kmer_hash_ungapped_test, issue1614)
{
    std::vector<seqan3::dna5> sequence{"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5};
    EXPECT_RANGE_EQ(sequence | seqan3::views::kmer_hash(seqan3::ungapped{25}),
                    seqan3::views::repeat_n(298'023'223'876'953'124, 26));
}

// https://github.com/seqan/seqan3/issues/1643
TEST(kmer_hash_ungapped_test, issue1643)
{
    std::vector<seqan3::dna4> text_23_elements{"ACGATCGATCGTAGCTACTGAGC"_dna4};

    auto k_mer_size_23_view = text_23_elements | seqan3::views::kmer_hash(seqan3::ungapped{23u});
    EXPECT_EQ(k_mer_size_23_view.size(), 1u);
    EXPECT_EQ(k_mer_size_23_view[0], 6'829'917'194'121u);

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

// https://github.com/seqan/seqan3/issues/1754
TYPED_TEST(kmer_hash_ungapped_test, issue1754)
{
    TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC

    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        EXPECT_RANGE_EQ(result_t{36}, text1 | prefix_until_first_thymine | std::views::reverse | ungapped_view);
    }
}

// https://github.com/seqan/seqan3/issues/1754
TYPED_TEST(kmer_hash_gapped_test, issue1754)
{
    TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC

    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        EXPECT_RANGE_EQ(result_t{8}, text1 | prefix_until_first_thymine | std::views::reverse | gapped_view);
    }
}

// https://github.com/seqan/seqan3/issues/1953
TEST(kmer_hash_ungapped_test, issue1953)
{
    std::vector<seqan3::dna4> text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
    auto v = text1 | prefix_until_first_thymine | std::views::reverse | ungapped_view;
    EXPECT_EQ(1u, v.size());
}

// https://github.com/seqan/seqan3/issues/1953
TEST(kmer_hash_gapped_test, issue1953)
{
    std::vector<seqan3::dna4> text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
    auto v = text1 | prefix_until_first_thymine | std::views::reverse | gapped_view;
    EXPECT_EQ(1u, v.size());
}

// https://github.com/seqan/seqan3/issues/1963
TYPED_TEST(kmer_hash_ungapped_test, issue1963)
{
    TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
    result_t ungapped{57, 36, 19, 13, 54};
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        auto v = text1 | seqan3::views::complement | ungapped_view;
        EXPECT_RANGE_EQ(ungapped, v);
        EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    }
}

// https://github.com/seqan/seqan3/issues/1963
TYPED_TEST(kmer_hash_gapped_test, issue1963)
{
    TypeParam text1{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
    result_t gapped{13, 8, 7, 1, 14};
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        auto v = text1 | seqan3::views::complement | gapped_view;
        EXPECT_RANGE_EQ(gapped, v);
        EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    }
}

// https://github.com/seqan/seqan3/issues/1988
TYPED_TEST(kmer_hash_ungapped_test, issue1988)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
        result_t ungapped{6, 27, 44, 50, 9};

        auto v = text | ungapped_view;
        EXPECT_RANGE_EQ(ungapped, v);
        EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);

        auto v2 = text | ungapped_view | std::views::reverse;
        EXPECT_RANGE_EQ(ungapped | std::views::reverse, v2);
        EXPECT_TRUE(seqan3::const_iterable_range<decltype(v2)>);
    }
}

// https://github.com/seqan/seqan3/issues/1988
TYPED_TEST(kmer_hash_gapped_test, issue1988)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4}; // ACGTAGC
        result_t gapped{2, 7, 8, 14, 1};

        auto v = text | gapped_view;
        EXPECT_RANGE_EQ(gapped, v);
        EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);

        auto v2 = text | gapped_view | std::views::reverse;
        EXPECT_RANGE_EQ(gapped | std::views::reverse, v2);
        EXPECT_TRUE(seqan3::const_iterable_range<decltype(v2)>);
    }
}

// https://github.com/seqan/seqan3/issues/2415
TYPED_TEST(kmer_hash_ungapped_test, issue2415)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        TypeParam text{'T'_dna4, 'A'_dna4, 'A'_dna4}; // TAA
        result_t ungapped{48};

        auto v = text | ungapped_view | std::views::reverse;
        EXPECT_RANGE_EQ(ungapped, v);
    }
}

// https://github.com/seqan/seqan3/issues/2415
TYPED_TEST(kmer_hash_gapped_test, issue2415)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        TypeParam text{'T'_dna4, 'A'_dna4, 'A'_dna4}; // TAA
        result_t gapped{12};

        auto v = text | gapped_view | std::views::reverse;
        EXPECT_RANGE_EQ(gapped, v);
    }
}
