// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms fi the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/alignment_file/detail.hpp>

using namespace seqan3;

template <typename T>
class aligned_sequence : public ::testing::Test
{};

TYPED_TEST_CASE_P(aligned_sequence);

TYPED_TEST_P(aligned_sequence, fulfills_concept)
{
    EXPECT_TRUE((aligned_sequence_concept<TypeParam>));
}

TYPED_TEST_P(aligned_sequence, insert_one_gap)
{
    TypeParam aligned_seq;
    TypeParam aligned_seq_expected;
    TestFixture::initialize_typed_test_container(aligned_seq, "ACTA");
    TestFixture::initialize_typed_test_container(aligned_seq_expected, "A-CTA");

    auto it = insert_gap(aligned_seq, begin(aligned_seq) + 1);

    EXPECT_EQ(*it, gap{});
    EXPECT_EQ(aligned_seq, aligned_seq_expected);
}

TYPED_TEST_P(aligned_sequence, insert_multiple_gaps)
{
    TypeParam aligned_seq;
    TypeParam aligned_seq_expected;
    TestFixture::initialize_typed_test_container(aligned_seq, "ACTA");
    TestFixture::initialize_typed_test_container(aligned_seq_expected, "A--CTA");

    auto it = insert_gap(aligned_seq, begin(aligned_seq) + 1, 2);

    EXPECT_EQ(*it, gap{});
    EXPECT_EQ(*++it, gap{});
    EXPECT_EQ(aligned_seq, aligned_seq_expected);
}

TYPED_TEST_P(aligned_sequence, erase_one_gap)
{
    // 1) Removing an actual gap
    TypeParam aligned_seq;
    TypeParam aligned_seq_expected;
    TestFixture::initialize_typed_test_container(aligned_seq, "A-CTA");
    TestFixture::initialize_typed_test_container(aligned_seq_expected, "ACTA");

    auto it = erase_gap(aligned_seq, begin(aligned_seq) + 1);

    typename TypeParam::value_type val{'C'_dna4};
    val = 'C'_dna4;
    EXPECT_EQ(*it, val);
    EXPECT_EQ(aligned_seq, aligned_seq_expected);

    // 2) Removing a non-gap
    TypeParam aligned_seq_fail;
    TestFixture::initialize_typed_test_container(aligned_seq_fail, "A-CTA");

    EXPECT_THROW(erase_gap(aligned_seq_fail, begin(aligned_seq_fail) + 2), gap_erase_failure);

    EXPECT_EQ(aligned_seq_fail[1], gap{}); // nothing has changed
}

TYPED_TEST_P(aligned_sequence, erase_multiple_gaps)
{
    // 1) Removing actual gaps

    // nucleotide alphabet
    TypeParam aligned_seq;
    TypeParam aligned_seq_expected;
    TestFixture::initialize_typed_test_container(aligned_seq, "A--CTA");
    TestFixture::initialize_typed_test_container(aligned_seq_expected, "ACTA");

    auto it = erase_gap(aligned_seq, begin(aligned_seq) + 1, begin(aligned_seq) + 3);

    typename TypeParam::value_type val{'C'_dna4};
    val = 'C'_dna4;
    EXPECT_EQ(*it, val);
    EXPECT_EQ(aligned_seq, aligned_seq_expected);

    // 2) Removing a non-gap
    TypeParam aligned_seq_fail;
    TestFixture::initialize_typed_test_container(aligned_seq_fail, "A-CTA");

    EXPECT_THROW(erase_gap(aligned_seq_fail, begin(aligned_seq_fail) + 1, begin(aligned_seq_fail) + 3),
                 gap_erase_failure);

    EXPECT_EQ(aligned_seq_fail[1], gap{}); // nothing has changed
}

TYPED_TEST_P(aligned_sequence, cigar_string)
{
    {   // default_parameters
        TypeParam ref;
        TypeParam read;
        TestFixture::initialize_typed_test_container(ref,  "ACGTGAT--CTG");
        TestFixture::initialize_typed_test_container(read, "ACGT-CGTAGTG");

        std::string expected = "4M1D2M2I3M";

        EXPECT_EQ(expected, detail::get_cigar_string(std::make_pair(ref, read)));

        TypeParam ref2;
        TypeParam read2;
        TestFixture::initialize_typed_test_container(ref2,  "---ACGTGAT--CTG--");
        TestFixture::initialize_typed_test_container(read2, "-ACGT-CGTAGTG----");

        std::string expected2 = "1P2I2M1D4M2I1M2D2P";

        EXPECT_EQ(expected2, detail::get_cigar_string(std::make_pair(ref2, read2)));
    }
    {
        // with soft clipping
        TypeParam ref;
        TypeParam read;
        TestFixture::initialize_typed_test_container(ref,  "ACGTGAT--CTG");
        TestFixture::initialize_typed_test_container(read, "ACGT-CGTAGTG");

        std::string expected = "5S4M1D2M2I3M60S";

        EXPECT_EQ(expected, detail::get_cigar_string(std::make_pair(ref, read), 5, 60));

        // gaps at the end
        TypeParam ref2;
        TypeParam read2;
        TestFixture::initialize_typed_test_container(ref2,  "---ACGTGAT--CTG--");
        TestFixture::initialize_typed_test_container(read2, "-ACGT-CGTAGTG----");

        std::string expected2 = "3S1P2I2M1D4M2I1M2D2P5S";

        EXPECT_EQ(expected2, detail::get_cigar_string(std::make_pair(ref2, read2), 3, 5));
    }
    {
        // no gaps at the end
        TypeParam ref;
        TypeParam read;
        TestFixture::initialize_typed_test_container(ref,  "ACGTGAT--CAG");
        TestFixture::initialize_typed_test_container(read, "ACGT-CGTACTG");

        std::string expected1 =    "4=1D2X2I1=1X1=";
        std::string expected2 = "5S4=1D2X2I1=1X1=60S";

        EXPECT_EQ(expected1, detail::get_cigar_string(std::make_pair(ref, read), 0, 0, true));
        EXPECT_EQ(expected2, detail::get_cigar_string(std::make_pair(ref, read), 5, 60, true));
    }
}

REGISTER_TYPED_TEST_CASE_P(aligned_sequence, fulfills_concept, insert_one_gap, insert_multiple_gaps,
                           erase_one_gap, erase_multiple_gaps, cigar_string);
