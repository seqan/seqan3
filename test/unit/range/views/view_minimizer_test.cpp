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
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimizer.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>

#include <gtest/gtest.h>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_shape;

class minimizer_test : public ::testing::Test
{
protected:
    using result_t = std::vector<size_t>;

    static constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
    static constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);
    static constexpr auto minimizer_view = seqan3::views::minimizer(5);
    static constexpr auto minimizer_view2 = seqan3::views::minimizer(1); // kmer_size == window_size

    std::vector<seqan3::dna4> text1{"AAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same result for ungapped and gapped

    std::vector<seqan3::dna4> text2{"ACGTCGACGTTTAG"_dna4};
    result_t ungapped_no_rev2{27, 97, 27};   // ACGT, CGAC, ACGT
    result_t gapped_no_rev2{3, 5, 3};        // A--T, C--C, A--T - "-" for gap

    seqan3::bitcompressed_vector<seqan3::dna4> text3{"AC"_dna4};
    result_t result3{}; // Same result for ungapped and gapped

    seqan3::bitcompressed_vector<seqan3::dna4> text4{"ACGGCGACGTTTAG"_dna4};
    result_t ungapped_no_rev4{26, 97, 27};    // ACGG, CGAC, ACGT
    result_t gapped_no_rev4{2, 5, 3};         // A--G, C--C-, A--T "-" for gap
    result_t ungapped_no_rev4_stop{26, 97};   // For stop at first T
    result_t gapped_no_rev4_stop{2, 5};
};

TEST_F(minimizer_test, ungapped)
{
    EXPECT_EQ(result1, text1 | kmer_view | minimizer_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result1, text1_short | kmer_view | minimizer_view2 | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev2, text2 | kmer_view | minimizer_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result3, text3 | seqan3::views::kmer_hash(seqan3::ungapped{3}) | seqan3::views::to<result_t>);
    EXPECT_EQ(ungapped_no_rev4, text4 | kmer_view | minimizer_view | seqan3::views::to<result_t>);

}

TEST_F(minimizer_test, gapped)
{
    EXPECT_EQ(result1, text1 | gapped_kmer_view | minimizer_view | seqan3::views::to<result_t>);
    EXPECT_EQ(result1, text1_short | gapped_kmer_view | minimizer_view2 | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev2, text2 | gapped_kmer_view | minimizer_view | seqan3::views::to<result_t>);
    auto test =  text3 | seqan3::views::kmer_hash(0b101_shape);
    seqan3::debug_stream << test;
    auto test2 = test | minimizer_view;
    seqan3::debug_stream << test2;
    EXPECT_EQ(result3, text3 | seqan3::views::kmer_hash(0b101_shape) | minimizer_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev4, text4 | gapped_kmer_view | minimizer_view | seqan3::views::to<result_t>);
}

TEST_F(minimizer_test, combinability)
{
    auto stop_at_t = seqan3::views::take_until([] (seqan3::dna4 const x) { return x == 'T'_dna4; });
    EXPECT_EQ(ungapped_no_rev4_stop, text4 | stop_at_t | kmer_view | minimizer_view | seqan3::views::to<result_t>);
    EXPECT_EQ(gapped_no_rev4_stop, text4 | stop_at_t | gapped_kmer_view | minimizer_view | seqan3::views::to<result_t>);
}
/*
TEST_F(minimizer_test, invalid_sizes)
{
    EXPECT_NO_THROW(text1 | seqan3::views::minimizer(seqan3::ungapped{32}));
    EXPECT_THROW(text1 | seqan3::views::minimizer(seqan3::ungapped{33}), std::invalid_argument);
    EXPECT_NO_THROW(text1 | std::views::reverse | seqan3::views::minimizer(seqan3::ungapped{32}));
    EXPECT_THROW(text1 | std::views::reverse | seqan3::views::minimizer(seqan3::ungapped{33}), std::invalid_argument);

    EXPECT_NO_THROW(text1 | seqan3::views::minimizer(0xFFFFFFFE001_shape)); // size=44, count=32
    EXPECT_THROW(text1 | seqan3::views::minimizer(0xFFFFFFFFE009_shape), std::invalid_argument); // size=44, count=33

    std::vector<seqan3::dna5> dna5_text{};
    EXPECT_NO_THROW(dna5_text | seqan3::views::minimizer(seqan3::ungapped{27}));
    EXPECT_THROW(dna5_text | seqan3::views::minimizer(seqan3::ungapped{28}), std::invalid_argument);
}*/
