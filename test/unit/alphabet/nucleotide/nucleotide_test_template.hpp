// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

using namespace seqan3;

template <typename T>
class nucleotide : public ::testing::Test
{};

TYPED_TEST_CASE_P(nucleotide);

TYPED_TEST_P(nucleotide, global_complement)
{
    EXPECT_EQ(complement(TypeParam::A), TypeParam::T);
    EXPECT_EQ(complement(TypeParam::C), TypeParam::G);
    EXPECT_EQ(complement(TypeParam::G), TypeParam::C);
    EXPECT_EQ(complement(TypeParam::T), TypeParam::A);

    using vsize_t = std::decay_t<decltype(alphabet_size_v<TypeParam>)>;

    for (vsize_t i = 0u; i < alphabet_size_v<TypeParam>; ++i)
    {
        TypeParam c = assign_rank(TypeParam{}, i);

        EXPECT_EQ(complement(complement(c)), c);
    }
}

TYPED_TEST_P(nucleotide, concept_check)
{
    EXPECT_TRUE(nucleotide_concept<TypeParam>);
    EXPECT_TRUE(nucleotide_concept<TypeParam &>);
}

REGISTER_TYPED_TEST_CASE_P(nucleotide, global_complement, concept_check);
