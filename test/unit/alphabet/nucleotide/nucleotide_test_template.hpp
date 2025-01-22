// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

template <typename t>
using nucleotide = ::testing::Test;

TYPED_TEST_SUITE_P(nucleotide);

TYPED_TEST_P(nucleotide, concept_check)
{
    EXPECT_TRUE(seqan3::trivial<TypeParam>);

    EXPECT_TRUE(seqan3::nucleotide_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<TypeParam &>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<TypeParam const &>);
}

TYPED_TEST_P(nucleotide, complement)
{
    EXPECT_EQ(seqan3::complement(TypeParam{}.assign_char('A')), TypeParam{}.assign_char('T'));
    EXPECT_EQ(seqan3::complement(TypeParam{}.assign_char('C')), TypeParam{}.assign_char('G'));
    EXPECT_EQ(seqan3::complement(TypeParam{}.assign_char('G')), TypeParam{}.assign_char('C'));
    EXPECT_EQ(seqan3::complement(TypeParam{}.assign_char('T')), TypeParam{}.assign_char('A'));

    using vsize_t = std::decay_t<decltype(seqan3::alphabet_size<TypeParam>)>;

    for (vsize_t i = 0u; i < seqan3::alphabet_size<TypeParam>; ++i)
    {
        TypeParam c = seqan3::assign_rank_to(i, TypeParam{});

        EXPECT_EQ(seqan3::complement(seqan3::complement(c)), c);
    }
}

REGISTER_TYPED_TEST_SUITE_P(nucleotide, concept_check, complement);
