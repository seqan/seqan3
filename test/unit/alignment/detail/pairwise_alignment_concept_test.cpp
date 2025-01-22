// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/detail/pairwise_alignment_concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

TEST(pairwise_alignment_concept, std_pair_gapped_sequences)
{
    using gapped_sequence1_t = std::vector<seqan3::gapped<seqan3::dna4>>;
    using gapped_sequence2_t = std::vector<seqan3::gapped<seqan3::dna4>>;
    using alignment_t = std::pair<gapped_sequence1_t, gapped_sequence2_t>;

    EXPECT_FALSE((seqan3::detail::pairwise_alignment<gapped_sequence1_t>));
    EXPECT_FALSE((seqan3::detail::pairwise_alignment<gapped_sequence2_t>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t &>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const &>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t &&>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const &&>));
}

TEST(pairwise_alignment_concept, std_tuple_gapped_sequences)
{
    using gapped_sequence1_t = std::vector<seqan3::gapped<seqan3::dna4>> const;
    using gapped_sequence2_t = std::vector<seqan3::gapped<seqan3::dna5>>;
    using alignment_t = std::tuple<gapped_sequence1_t, gapped_sequence2_t>;

    EXPECT_FALSE((seqan3::detail::pairwise_alignment<gapped_sequence1_t>));
    EXPECT_FALSE((seqan3::detail::pairwise_alignment<gapped_sequence2_t>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t &>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const &>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t &&>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const &&>));
}

TEST(pairwise_alignment_concept, std_tuple_gap_sequence)
{
    using gap_sequence1_t = std::vector<seqan3::gap> const;
    using gap_sequence2_t = std::vector<seqan3::gap>;
    using alignment_t = std::tuple<gap_sequence1_t, gap_sequence2_t>;

    EXPECT_FALSE((seqan3::detail::pairwise_alignment<gap_sequence1_t>));
    EXPECT_FALSE((seqan3::detail::pairwise_alignment<gap_sequence2_t>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t &>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const &>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t &&>));
    EXPECT_TRUE((seqan3::detail::pairwise_alignment<alignment_t const &&>));
}

TEST(writable_pairwise_alignment_concept, std_pair_gapped_sequences)
{
    using gapped_sequence1_t = std::vector<seqan3::gapped<seqan3::dna4>>;
    using gapped_sequence2_t = std::vector<seqan3::gapped<seqan3::dna4>>;
    using alignment_t = std::pair<gapped_sequence1_t, gapped_sequence2_t>;

    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<gapped_sequence1_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<gapped_sequence2_t>));
    EXPECT_TRUE((seqan3::detail::writable_pairwise_alignment<alignment_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const>));
    EXPECT_TRUE((seqan3::detail::writable_pairwise_alignment<alignment_t &>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const &>));
    EXPECT_TRUE((seqan3::detail::writable_pairwise_alignment<alignment_t &&>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const &&>));
}

TEST(writable_pairwise_alignment_concept, std_tuple_gapped_sequences)
{
    using gapped_sequence1_t = std::vector<seqan3::gapped<seqan3::dna4>>;
    using gapped_sequence2_t = std::vector<seqan3::gapped<seqan3::dna5>> const;
    using alignment_t = std::tuple<gapped_sequence1_t, gapped_sequence2_t>;

    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<gapped_sequence1_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<gapped_sequence2_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t &>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const &>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t &&>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const &&>));
}

TEST(writable_pairwise_alignment_concept, std_tuple_gap_sequence)
{
    using gap_sequence1_t = std::vector<seqan3::gap>;
    using gap_sequence2_t = std::vector<seqan3::gap> const;
    using alignment_t = std::tuple<gap_sequence1_t, gap_sequence2_t>;

    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<gap_sequence1_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<gap_sequence2_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t &>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const &>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t &&>));
    EXPECT_FALSE((seqan3::detail::writable_pairwise_alignment<alignment_t const &&>));
}
