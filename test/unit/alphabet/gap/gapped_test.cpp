// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

using gapped_types = ::testing::Types<gapped<dna4>, gapped<dna15>, gapped<qualified<dna4, phred42>>>;

INSTANTIATE_TYPED_TEST_CASE_P(gapped, alphabet, gapped_types);
INSTANTIATE_TYPED_TEST_CASE_P(gapped, alphabet_constexpr, gapped_types);

template <typename t>
using gapped_test = ::testing::Test;

TYPED_TEST_CASE(gapped_test, gapped_types);

TYPED_TEST(gapped_test, concept_check)
{
    EXPECT_TRUE((AlignedSequence<std::vector<TypeParam>>));
}

TEST(gapped_test, initialise_from_component_alphabet)
{
    using alphabet_t = gapped<dna4>;

    constexpr alphabet_t letter0{'A'_dna4};
    constexpr alphabet_t letter1 = 'C'_dna4;
    constexpr alphabet_t letter2 = {'G'_dna4};
    constexpr alphabet_t letter3 = static_cast<alphabet_t>('T'_dna4);

    alphabet_t letter4{'A'_dna4};
    alphabet_t letter5 = 'C'_dna4;
    alphabet_t letter6 = {'G'_dna4};
    alphabet_t letter7 = static_cast<alphabet_t>('T'_dna4);

    constexpr alphabet_t letter8{gap{}}; // letter3 = 'T'_dna4; does not work
    alphabet_t letter9{gap{}};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 0);
    EXPECT_EQ(letter5.to_rank(), 1);
    EXPECT_EQ(letter6.to_rank(), 2);
    EXPECT_EQ(letter7.to_rank(), 3);
    EXPECT_EQ(letter8.to_rank(), 4);
    EXPECT_EQ(letter9.to_rank(), 4);
}

TEST(gapped_test, assign_from_component_alphabet)
{
    using alphabet_t = gapped<dna4>;
    alphabet_t letter{};

    letter = gap{};
    EXPECT_EQ(letter.to_rank(), 4);

    letter = 'A'_dna4;
    EXPECT_EQ(letter.to_rank(), 0);

    letter = {'C'_dna4}; // letter = {'C'_dna4}; does not work
    EXPECT_EQ(letter.to_rank(), 1);

    letter = static_cast<alphabet_t>('G'_dna4);
    EXPECT_EQ(letter.to_rank(), 2);

    letter = {static_cast<alphabet_t>('T'_dna4)};
    EXPECT_EQ(letter.to_rank(), 3);
}
