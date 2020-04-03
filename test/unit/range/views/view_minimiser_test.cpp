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
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser.hpp>
#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>

#include <gtest/gtest.h>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_shape;

class minimiser_test : public ::testing::Test
{
protected:
    using result_t = std::vector<size_t>;

    static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
    static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);
    static constexpr auto minimiser_view = seqan3::views::minimiser(5);
    static constexpr auto minimiser_view2 = seqan3::views::minimiser(1); // kmer_size == window_size

    std::vector<seqan3::dna4> text1{"AAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text2{"ACGTCGACGTTTAG"_dna4};
    std::vector<seqan3::dna4> ctext2{"ACGTCGACGTTTAG"_dna4};
    seqan3::bitcompressed_vector<seqan3::dna4> bit_text2{"ACGTCGACGTTTAG"_dna4};
    seqan3::bitcompressed_vector<seqan3::dna4> const bit_ctext2{"ACGTCGACGTTTAG"_dna4};
    std::list<seqan3::dna4> list_text2{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4, 'C'_dna4,
                                        'G'_dna4, 'T'_dna4, 'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4};
    std::list<seqan3::dna4> const list_ctext2{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4,
                                              'C'_dna4, 'G'_dna4, 'T'_dna4, 'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4};
    std::forward_list<seqan3::dna4>  flist_text2{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4,
                                                 'C'_dna4, 'G'_dna4, 'T'_dna4, 'T'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4};
    std::forward_list<seqan3::dna4> const flist_ctext2{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'C'_dna4, 'G'_dna4,
                                                       'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'T'_dna4, 'T'_dna4,
                                                       'A'_dna4, 'G'_dna4};
    result_t ungapped_no_rev2{27, 97, 27};   // ACGT, CGAC, ACGT
    result_t gapped_no_rev2{3, 5, 3};        // A--T, C--C, A--T - "-" for gap


    std::vector<seqan3::dna4> text3{"AC"_dna4};
    result_t result3{}; // Same result for ungapped and gapped

    seqan3::bitcompressed_vector<seqan3::dna4> text4{"ACGGCGACGTTTAG"_dna4};
    result_t ungapped_no_rev4{26, 97, 27};    // ACGG, CGAC, ACGT
    result_t gapped_no_rev4{2, 5, 3};         // A--G, C--C-, A--T "-" for gap
    result_t ungapped_no_rev4_stop{26, 97};   // For stop at first T
    result_t gapped_no_rev4_stop{2, 5};
};

TEST_F(minimiser_test, concepts)
{
    auto v1 = text1 | kmer_view | minimiser_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));

    auto v2 = list_text2 | kmer_view | minimiser_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), size_t>));

    auto v3 = flist_text2 | kmer_view | minimiser_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::view<decltype(v3)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v3)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v3)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v3), size_t>));
}

TEST_F(minimiser_test, different_inputs_kmer_hash)
{
    EXPECT_EQ(ungapped_no_rev2, text2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, ctext2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, bit_text2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, bit_ctext2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, list_text2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, list_ctext2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, flist_text2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, flist_ctext2 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
}

TEST_F(minimiser_test, ungapped_kmer_hash)
{
    EXPECT_EQ(result1, text1 | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result1, text1_short | kmer_view | minimiser_view2 | seqan3::views::to<result_t>);
    // Can be changed to kmer_view once #1699 is merged
    EXPECT_EQ(result3, text3 | seqan3::views::kmer_hash(seqan3::ungapped{3}) | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev4, text4 | kmer_view | minimiser_view | seqan3::views::to<result_t>);

}

TEST_F(minimiser_test, gapped_kmer_hash)
{
    EXPECT_EQ(result1, text1 | gapped_kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result1, text1_short | gapped_kmer_view | minimiser_view2 | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev2, text2 | gapped_kmer_view | minimiser_view | seqan3::views::to<result_t>);
    // Can be changed to gapped_kmer_view once #1699 is merged
    EXPECT_EQ(result3, text3 | seqan3::views::kmer_hash(0b101_shape) | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev4, text4 | gapped_kmer_view | minimiser_view | seqan3::views::to<result_t>);
}

TEST_F(minimiser_test, combinability)
{
    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_EQ(ungapped_no_rev4_stop, text4 | stop_at_t | kmer_view | minimiser_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev4_stop, text4 | stop_at_t | gapped_kmer_view | minimiser_view | seqan3::views::to<result_t>);
}
