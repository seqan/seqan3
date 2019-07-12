// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <range/v3/view/iota.hpp>

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/char_operations/predicate.hpp>

using namespace seqan3;

template <typename T>
using nucleotide = ::testing::Test;

TYPED_TEST_CASE_P(nucleotide);

TYPED_TEST_P(nucleotide, concept_check)
{
    EXPECT_TRUE(NucleotideAlphabet<TypeParam>);
    EXPECT_TRUE(NucleotideAlphabet<TypeParam &>);
    EXPECT_TRUE(NucleotideAlphabet<TypeParam const>);
    EXPECT_TRUE(NucleotideAlphabet<TypeParam const &>);
}

TYPED_TEST_P(nucleotide, global_complement)
{
    EXPECT_EQ(complement(TypeParam{}.assign_char('A')), TypeParam{}.assign_char('T'));
    EXPECT_EQ(complement(TypeParam{}.assign_char('C')), TypeParam{}.assign_char('G'));
    EXPECT_EQ(complement(TypeParam{}.assign_char('G')), TypeParam{}.assign_char('C'));
    EXPECT_EQ(complement(TypeParam{}.assign_char('T')), TypeParam{}.assign_char('A'));

    using vsize_t = std::decay_t<decltype(alphabet_size<TypeParam>)>;

    for (vsize_t i = 0u; i < alphabet_size<TypeParam>; ++i)
    {
        TypeParam c = assign_rank_to(i, TypeParam{});

        EXPECT_EQ(complement(complement(c)), c);
    }
}

REGISTER_TYPED_TEST_CASE_P(nucleotide, concept_check, global_complement);
