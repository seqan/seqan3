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
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>

#include <gtest/gtest.h>

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{4};
static constexpr seqan3::shape gapped_shape = 0b1001_shape;
static constexpr auto ungapped_view = seqan3::views::minimiser_hash(ungapped_shape, 8, 0); // window_size 8, seed 0
static constexpr auto gapped_view = seqan3::views::minimiser_hash(gapped_shape, 8, 0);     // window_size 8, seed 0

template <typename T>
class minimiser_hash_properties_test: public ::testing::Test { };

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitcompressed_vector<seqan3::dna4>,
                                                seqan3::bitcompressed_vector<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;

TYPED_TEST_SUITE(minimiser_hash_properties_test, underlying_range_types, );

class minimiser_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same for ungapped and gapped
    result_t result1_default_seed{0x8F3F73B5CF1C9ADE, 0x8F3F73B5CF1C9ADE, 0x8F3F73B5CF1C9ADE}; // Default seed used

    std::vector<seqan3::dna4> text2{"AC"_dna4};
    result_t result2{};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t ungapped_no_rev3{26, 97}; // ACGG, CGAC
    result_t gapped_no_rev3{2, 5};     // A--G, C--C- "-" for gap
};

TYPED_TEST(minimiser_hash_properties_test, different_input_ranges)
{
    TypeParam text{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4,
                   'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4}; // ACGTCGACGTTTAG
    result_t ungapped_no_rev{27, 97, 27}; // ACGT, CGAC, ACGT
    result_t gapped_no_rev{3, 5, 3};      // A--T, C--C, A--T - "-" for gap
    EXPECT_EQ(ungapped_no_rev, text | ungapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev, text | gapped_view | seqan3::views::to<result_t>);
}

TEST_F(minimiser_hash_test, ungapped)
{
    EXPECT_EQ(result1, text1 | ungapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result2, text2 | ungapped_view | seqan3::views::to<result_t>);

    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_EQ(ungapped_no_rev3, text3 | stop_at_t | seqan3::views::minimiser_hash(seqan3::ungapped{4}, 8, 0)
                                      | seqan3::views::to<result_t>);
}

TEST_F(minimiser_hash_test, gapped)
{
    EXPECT_EQ(result1, text1 | gapped_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result2, text2 | gapped_view | seqan3::views::to<result_t>);

    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_EQ(gapped_no_rev3, text3 | stop_at_t | gapped_view | seqan3::views::to<result_t>);
}

TEST_F(minimiser_hash_test, default_seed)
{
    EXPECT_EQ(result1_default_seed, text1 | seqan3::views::minimiser_hash(ungapped_shape, 8)
                                          | seqan3::views::to<result_t>);
    EXPECT_EQ(result1_default_seed, text1 | seqan3::views::minimiser_hash(gapped_shape, 8)
                                          | seqan3::views::to<result_t>);
}

TEST_F(minimiser_hash_test, shape_size_equal_window_size)
{
    auto apply_seed = std::views::transform([] (int i) { return i ^ 0x8F3F73B5CF1C9ADE; });
    EXPECT_EQ(text1 | seqan3::views::kmer_hash(ungapped_shape)
                    | apply_seed
                    | seqan3::views::to<result_t>,
              text1 | seqan3::views::minimiser_hash(ungapped_shape) | seqan3::views::to<result_t>);
    EXPECT_EQ(text1 | seqan3::views::kmer_hash(gapped_shape)
                    | apply_seed
                    | seqan3::views::to<result_t>,
              text1 | seqan3::views::minimiser_hash(gapped_shape) | seqan3::views::to<result_t>);
}


TEST_F(minimiser_hash_test, shape_bigger_than_window)
{
    EXPECT_THROW(text1 | seqan3::views::minimiser_hash(ungapped_shape, 3), std::invalid_argument);
    EXPECT_THROW(text1 | seqan3::views::minimiser_hash(gapped_shape, 3), std::invalid_argument);
}
