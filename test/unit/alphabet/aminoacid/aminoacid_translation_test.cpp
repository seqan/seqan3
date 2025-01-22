// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using seqan3::operator""_aa27;
using seqan3::operator""_dna4;
using seqan3::operator""_dna15;

TEST(translate_triplets, dna4)
{
    seqan3::dna4 n1{'C'_dna4};
    seqan3::dna4 n2{'T'_dna4};
    seqan3::dna4 n3{'A'_dna4};
    seqan3::aa27 c{'L'_aa27};

    seqan3::aa27 t1{seqan3::translate_triplet<seqan3::genetic_code::canonical>(n1, n2, n3)};

    EXPECT_EQ(t1, c);
}

TEST(translate_triplets, dna15)
{
    seqan3::dna15 n1{'C'_dna15};
    seqan3::dna15 n2{'T'_dna15};
    seqan3::dna15 n3{'A'_dna15};
    seqan3::aa27 c{'L'_aa27};

    seqan3::aa27 t1{seqan3::translate_triplet<seqan3::genetic_code::canonical, seqan3::dna15>(n1, n2, n3)};

    EXPECT_EQ(t1, c);
}
