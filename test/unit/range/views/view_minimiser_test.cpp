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
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "../iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;
using iterator_type = std::ranges::iterator_t< decltype(std::declval<seqan3::dna4_vector&>()
                                                        | seqan3::views::kmer_hash(seqan3::ungapped{4})
                                                        | seqan3::views::minimiser(5))>;

static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);
static constexpr auto minimiser_view = seqan3::views::minimiser(5);
static constexpr auto minimiser_view2 = seqan3::views::minimiser(1); // kmer_size == window_size, should throw

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{26, 97, 27};

    decltype(seqan3::views::minimiser(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 5)) test_range =
    seqan3::views::minimiser(vec, 5);
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class minimiser_view_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitcompressed_vector<seqan3::dna4>,
                                                seqan3::bitcompressed_vector<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;

TYPED_TEST_SUITE(minimiser_view_properties_test, underlying_range_types, );

class minimiser_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same result for ungapped and gapped
    result_t result1_short{0}; // windows_size == text_size, same result for ungapped and gapped

    std::vector<seqan3::dna4> too_short_text{"AC"_dna4};

    seqan3::bitcompressed_vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t result3_ungapped_no_rev{26, 97, 27}; // ACGG, CGAC, ACGT
    result_t result3_gapped_no_rev{2, 5, 3};      // A--G, C--C-, A--T "-" for gap
    result_t result3_ungapped_no_rev_stop{26, 97};// For stop at first T
    result_t result3_gapped_no_rev_stop{2, 5};    // For stop at first T
};

TYPED_TEST(minimiser_view_properties_test, concepts)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    auto v = text | kmer_view | minimiser_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(minimiser_view_properties_test, different_inputs_kmer_hash)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped_no_rev{27, 97, 27};                // ACGT, CGAC, ACGT
    result_t gapped_no_rev{3, 5, 3};                     // A--T, C--C, A--T - "-" for gap
    EXPECT_RANGE_EQ(ungapped_no_rev, text | kmer_view | minimiser_view);
    EXPECT_RANGE_EQ(gapped_no_rev, text | gapped_kmer_view | minimiser_view);
}

TEST_F(minimiser_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | minimiser_view);
    EXPECT_THROW(text1_short | kmer_view | minimiser_view2, std::invalid_argument);
    auto empty_view = too_short_text | kmer_view | minimiser_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3_ungapped_no_rev, text3 | kmer_view | minimiser_view);

}

TEST_F(minimiser_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | minimiser_view);
    EXPECT_THROW(text1_short | gapped_kmer_view | minimiser_view2, std::invalid_argument);
    auto empty_view = too_short_text | gapped_kmer_view | minimiser_view;
    EXPECT_TRUE(std::ranges::empty(empty_view));
    EXPECT_RANGE_EQ(result3_gapped_no_rev, text3 | gapped_kmer_view | minimiser_view);
}

TEST_F(minimiser_test, window_too_big)
{
    EXPECT_RANGE_EQ(result1_short, text1 | kmer_view | seqan3::views::minimiser(20));
    EXPECT_RANGE_EQ(result1_short, text1 | gapped_kmer_view | seqan3::views::minimiser(20));
}

TEST_F(minimiser_test, combinability)
{
    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_RANGE_EQ(result3_ungapped_no_rev_stop, text3 | stop_at_t | kmer_view | minimiser_view);
    EXPECT_RANGE_EQ(result3_gapped_no_rev_stop, text3 | stop_at_t | gapped_kmer_view | minimiser_view);
}
