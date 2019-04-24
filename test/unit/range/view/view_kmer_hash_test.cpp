// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/view/kmer_hash.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

TEST(kmer_hash, view)
{
    {
        std::vector<dna4> text{"AAAAA"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{0,0,0};
        EXPECT_EQ(expected, hashes);
    }
    {
        std::vector<dna4> text{"ACGTAGC"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{6,27,44,50,9};
        EXPECT_EQ(expected, hashes);
    }
    {
        std::vector<dna4> text{"AC"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{};
        EXPECT_EQ(expected, hashes);
    }
    {
        bitcompressed_vector<dna4> text{"ACGTAGC"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{6,27,44,50,9};
        EXPECT_EQ(expected, hashes);
    }
}

TEST(kmer_hash, const_view)
{
    {
        std::vector<dna4> const text{"AAAAA"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{0,0,0};
        EXPECT_EQ(expected, hashes);
    }
    {
        std::vector<dna4> const text{"ACGTAGC"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{6,27,44,50,9};
        EXPECT_EQ(expected, hashes);
    }
    {
        std::vector<dna4> const text{"AC"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{};
        EXPECT_EQ(expected, hashes);
    }
    {
        bitcompressed_vector<dna4> const text{"ACGTAGC"_dna4};
        std::vector<size_t> hashes = text | view::kmer_hash(3);
        std::vector<size_t> expected{6,27,44,50,9};
        EXPECT_EQ(expected, hashes);
    }
}
