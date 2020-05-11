// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/concept.hpp>

template <typename T>
using phred = ::testing::Test;

TYPED_TEST_SUITE_P(phred);

// test provision of data type `phred_type` and phred converter.
TYPED_TEST_P(phred, concept_check)
{
    EXPECT_TRUE(seqan3::quality_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::quality_alphabet<TypeParam &>);
    EXPECT_TRUE(seqan3::quality_alphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::quality_alphabet<TypeParam const &>);

    EXPECT_TRUE(seqan3::writable_quality_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::writable_quality_alphabet<TypeParam &>);
    EXPECT_FALSE(seqan3::writable_quality_alphabet<TypeParam const>);
    EXPECT_FALSE(seqan3::writable_quality_alphabet<TypeParam const &>);
}

// more elaborate test of assign_char and to_char, basic test is in alphabet_test.cpp
TYPED_TEST_P(phred, conversion_char)
{
    using c_t = typename TypeParam::char_type;
    for (c_t i = std::numeric_limits<c_t>::lowest(); i < std::numeric_limits<c_t>::max(); ++i)
    {
        TypeParam v;
        v.assign_char(i);

        if (i < TypeParam::offset_char)                                     // too small, map to valid smallest
            EXPECT_EQ(v.to_char(), TypeParam::offset_char);
        else if (i >= TypeParam::offset_char + TypeParam::alphabet_size)    // too big, map to valid biggest
            EXPECT_EQ(v.to_char(), TypeParam::offset_char + TypeParam::alphabet_size - 1);
        else                                                                // valid range, map to identity
            EXPECT_EQ(v.to_char(), i);
    }
}

// test assign_phred and to_phred
TYPED_TEST_P(phred, conversion_phred)
{
    using p_t = typename TypeParam::phred_type;
    for (p_t i = std::numeric_limits<p_t>::lowest(); i < std::numeric_limits<p_t>::max(); ++i)
    {
        TypeParam v;
        v.assign_phred(i);

        if (i < TypeParam::offset_phred)                                    // too small, map to valid smallest
            EXPECT_EQ(v.to_phred(), TypeParam::offset_phred);
        else if (i >= TypeParam::offset_phred + TypeParam::alphabet_size)   // too big, map to valid biggest
            EXPECT_EQ(v.to_phred(), TypeParam::offset_phred + TypeParam::alphabet_size - 1);
        else                                                                // valid range, map to identity
            EXPECT_EQ(v.to_phred(), i);
    }
}

// test user-defined constructor
TYPED_TEST_P(phred, conversion_rank)
{
    TypeParam v{0};
    EXPECT_EQ(v.to_phred(), 0);
    EXPECT_EQ(v.to_rank(),  -TypeParam::offset_phred);

    TypeParam v2{23};
    EXPECT_EQ(v2.to_phred(), 23);
    EXPECT_EQ(v2.to_rank(),  23 - TypeParam::offset_phred);
}

REGISTER_TYPED_TEST_SUITE_P(phred, concept_check, conversion_char, conversion_phred, conversion_rank);
