// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/scoring/hamming_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

TEST(hamming_scoring_scheme_test, default_constructor)
{
    seqan3::hamming_scoring_scheme scheme{};
    EXPECT_EQ(scheme.score('A', 'A'), 0);
    EXPECT_EQ(scheme.score('A', 'C'), -1);
    EXPECT_EQ(scheme.score('A', 'G'), -1);
    EXPECT_EQ(scheme.score('A', 'T'), -1);
}

TEST(hamming_scoring_scheme_test, models_scoring_scheme)
{
    EXPECT_TRUE((seqan3::scoring_scheme_for<seqan3::hamming_scoring_scheme, seqan3::dna4, seqan3::dna4>));
    EXPECT_TRUE((seqan3::scoring_scheme_for<seqan3::hamming_scoring_scheme, char, char>));
    EXPECT_TRUE((seqan3::scoring_scheme_for<seqan3::hamming_scoring_scheme, uint8_t, char>));
    EXPECT_TRUE((seqan3::scoring_scheme_for<seqan3::hamming_scoring_scheme, uint8_t, uint8_t>));
}

TEST(hamming_scoring_scheme_test, comparison)
{
    seqan3::hamming_scoring_scheme scheme{};
    EXPECT_TRUE(scheme == seqan3::hamming_scoring_scheme{});
    EXPECT_FALSE(scheme != seqan3::hamming_scoring_scheme{});
}
