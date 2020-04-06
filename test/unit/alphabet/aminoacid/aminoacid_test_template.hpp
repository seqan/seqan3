// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/core/char_operations/predicate.hpp>

template <typename T>
using aminoacid = ::testing::Test;

TYPED_TEST_SUITE_P(aminoacid);

TYPED_TEST_P(aminoacid, concept_check)
{
    EXPECT_TRUE(seqan3::aminoacid_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::aminoacid_alphabet<TypeParam &>);
    EXPECT_TRUE(seqan3::aminoacid_alphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::aminoacid_alphabet<TypeParam const &>);
}

// ------------------------------------------------------------------
// comparators
// ------------------------------------------------------------------

TYPED_TEST_P(aminoacid, comparators)
{
    EXPECT_TRUE(TypeParam{}.assign_char('A') == TypeParam{}.assign_char('A'));
    EXPECT_TRUE(TypeParam{}.assign_char('A') != TypeParam{}.assign_char('B'));
    EXPECT_TRUE(TypeParam{}.assign_char('A') <  TypeParam{}.assign_char('B'));
    EXPECT_TRUE(TypeParam{}.assign_char('A') <= TypeParam{}.assign_char('B'));
    EXPECT_TRUE(TypeParam{}.assign_char('B') >  TypeParam{}.assign_char('A'));
    EXPECT_TRUE(TypeParam{}.assign_char('B') >= TypeParam{}.assign_char('A'));
}

REGISTER_TYPED_TEST_SUITE_P(aminoacid, concept_check, comparators);
