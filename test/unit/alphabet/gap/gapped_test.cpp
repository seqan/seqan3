// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

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

TEST(gapped_test, initialise_from_component_alphabet)
{
    using alphabet_t = gapped<dna4>;

    constexpr alphabet_t letter0{dna4::A};
    constexpr alphabet_t letter1 = dna4::C;
    constexpr alphabet_t letter2 = {dna4::G};
    constexpr alphabet_t letter3 = static_cast<alphabet_t>(dna4::T);

    alphabet_t letter4{dna4::A};
    alphabet_t letter5 = dna4::C;
    alphabet_t letter6 = {dna4::G};
    alphabet_t letter7 = static_cast<alphabet_t>(dna4::T);

    constexpr alphabet_t letter8{gap::GAP}; // letter3 = dna4::T; does not work
    alphabet_t letter9{gap::GAP};

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

    letter = dna4::A;
    EXPECT_EQ(letter.to_rank(), 0);

    letter = {dna4::C}; // letter = {dna4::C}; does not work
    EXPECT_EQ(letter.to_rank(), 1);

    letter = static_cast<alphabet_t>(dna4::G);
    EXPECT_EQ(letter.to_rank(), 2);

    letter = {static_cast<alphabet_t>(dna4::T)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = gap::GAP;
    EXPECT_EQ(letter.to_rank(), 4);
}

TEST(gapped_test, fulfills_concepts)
{
    using alphabet_t = gapped<dna4>;
    EXPECT_TRUE((std::is_pod_v<alphabet_t>));
    EXPECT_TRUE((std::is_trivial_v<alphabet_t>));
    EXPECT_TRUE((std::is_trivially_copyable_v<alphabet_t>));
    EXPECT_TRUE((std::is_standard_layout_v<alphabet_t>));
}
