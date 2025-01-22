// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/cigar_conversion/alignment_from_cigar.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

struct cigar_conversion_data : public ::testing::Test
{
    seqan3::dna5_vector ref = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::vector<seqan3::cigar> simple_cigar{{1, 'S'_cigar_operation}, // 1S1M1D1M1I
                                            {1, 'M'_cigar_operation},
                                            {1, 'D'_cigar_operation},
                                            {1, 'M'_cigar_operation},
                                            {1, 'I'_cigar_operation}};

    std::vector<seqan3::gapped<seqan3::dna5>> simple_cigar_gapped_ref{'A'_dna5, 'C'_dna5, 'T'_dna5, seqan3::gap{}};
    std::vector<seqan3::gapped<seqan3::dna5>> simple_cigar_gapped_seq{'C'_dna5, seqan3::gap{}, 'G'_dna5, 'T'_dna5};

    std::vector<seqan3::cigar> cigar_with_padding{{1, 'S'_cigar_operation}, // 1S1M1P1M1I1M1I1D1M1S
                                                  {1, 'M'_cigar_operation},
                                                  {1, 'P'_cigar_operation},
                                                  {1, 'M'_cigar_operation},
                                                  {1, 'I'_cigar_operation},
                                                  {1, 'M'_cigar_operation},
                                                  {1, 'I'_cigar_operation},
                                                  {1, 'D'_cigar_operation},
                                                  {1, 'M'_cigar_operation},
                                                  {1, 'S'_cigar_operation}};

    // same as `cigar_with_padding` but M are substituted by `=` or `X` depending on match or mismatch
    std::vector<seqan3::cigar> extended_cigar_with_padding{{1, 'S'_cigar_operation}, // 1S1=1P1X1I1X1I1D1=1S
                                                           {1, '='_cigar_operation},
                                                           {1, 'P'_cigar_operation},
                                                           {1, 'X'_cigar_operation},
                                                           {1, 'I'_cigar_operation},
                                                           {1, 'X'_cigar_operation},
                                                           {1, 'I'_cigar_operation},
                                                           {1, 'D'_cigar_operation},
                                                           {1, '='_cigar_operation},
                                                           {1, 'S'_cigar_operation}};

    std::vector<seqan3::gapped<seqan3::dna5>> cigar_with_padding_gapped_ref =
        {'T'_dna5, seqan3::gap{}, 'G'_dna5, seqan3::gap{}, 'A'_dna5, seqan3::gap{}, 'T'_dna5, 'C'_dna5};
    std::vector<seqan3::gapped<seqan3::dna5>> cigar_with_padding_gapped_seq =
        {'T'_dna5, seqan3::gap{}, 'A'_dna5, 'G'_dna5, 'T'_dna5, 'A'_dna5, seqan3::gap{}, 'C'_dna5};

    std::vector<seqan3::cigar> cigar_with_hard_clipping{{1, 'H'_cigar_operation}, // 1H2S7M1D1M1S2H
                                                        {2, 'S'_cigar_operation},
                                                        {7, 'M'_cigar_operation},
                                                        {1, 'D'_cigar_operation},
                                                        {1, 'M'_cigar_operation},
                                                        {1, 'S'_cigar_operation},
                                                        {2, 'H'_cigar_operation}};

    std::vector<seqan3::gapped<seqan3::dna5>> cigar_with_hard_clipping_gapped_ref =
        {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5, 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5};
    std::vector<seqan3::gapped<seqan3::dna5>> cigar_with_hard_clipping_gapped_seq =
        {'A'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5, 'G'_dna5, 'N'_dna5, seqan3::gap{}, 'A'_dna5};
};

struct alignment_from_cigar : public cigar_conversion_data
{};

TEST_F(alignment_from_cigar, empty_cigar)
{
    std::vector<seqan3::cigar> empty_cigar{};
    seqan3::dna5_vector seq{"ACGT"_dna5};

    // An empty CIGAR string is not valid as it must always fulfil the following:
    // "Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ"
    EXPECT_THROW(seqan3::alignment_from_cigar(empty_cigar, this->ref, 0, seq), std::logic_error);
}

TEST_F(alignment_from_cigar, cigar_covers_too_little_bases_of_read)
{
    seqan3::dna5_vector seq{"ACGT"_dna5};
    std::vector<seqan3::cigar> corrupt_cigar{{3, 'M'_cigar_operation}}; // Although seq is of length 4

    EXPECT_THROW(seqan3::alignment_from_cigar(corrupt_cigar, this->ref, 0, seq), std::logic_error);
}

TEST_F(alignment_from_cigar, cigar_covers_too_many_bases_of_read)
{
    seqan3::dna5_vector seq{"ACGT"_dna5};
    std::vector<seqan3::cigar> corrupt_cigar{{5, 'M'_cigar_operation}}; // Although seq is of length 4

    EXPECT_THROW(seqan3::alignment_from_cigar(corrupt_cigar, this->ref, 0, seq), std::logic_error);
}

TEST_F(alignment_from_cigar, cigar_covers_too_many_bases_of_reference)
{
    seqan3::dna5_vector seq{"ACGT"_dna5};
    std::vector<seqan3::cigar> corrupt_cigar{{2, 'M'_cigar_operation},
                                             {40, 'D'_cigar_operation}, // Although reference is only of length 35
                                             {2, 'M'_cigar_operation}};

    EXPECT_THROW(seqan3::alignment_from_cigar(corrupt_cigar, this->ref, 0, seq), std::logic_error);
}

TEST_F(alignment_from_cigar, simple_cigar)
{
    seqan3::dna5_vector seq{"ACGT"_dna5};
    int32_t reference_start_position{0};

    auto alignment = seqan3::alignment_from_cigar(this->simple_cigar, this->ref, reference_start_position, seq);

    EXPECT_RANGE_EQ(std::get<0>(alignment), this->simple_cigar_gapped_ref);
    EXPECT_RANGE_EQ(std::get<1>(alignment), this->simple_cigar_gapped_seq);
}

TEST_F(alignment_from_cigar, with_padding)
{
    seqan3::dna5_vector seq{"GTAGTACA"_dna5};
    int32_t reference_start_position{2};

    auto alignment = seqan3::alignment_from_cigar(this->cigar_with_padding, this->ref, reference_start_position, seq);

    EXPECT_RANGE_EQ(std::get<0>(alignment), this->cigar_with_padding_gapped_ref);
    EXPECT_RANGE_EQ(std::get<1>(alignment), this->cigar_with_padding_gapped_seq);
}

TEST_F(alignment_from_cigar, extended_cigar)
{
    seqan3::dna5_vector seq{"GTAGTACA"_dna5};
    int32_t reference_start_position{2};

    auto alignment =
        seqan3::alignment_from_cigar(this->extended_cigar_with_padding, this->ref, reference_start_position, seq);

    EXPECT_RANGE_EQ(std::get<0>(alignment), this->cigar_with_padding_gapped_ref);
    EXPECT_RANGE_EQ(std::get<1>(alignment), this->cigar_with_padding_gapped_seq);
}

TEST_F(alignment_from_cigar, with_hard_clipping)
{
    seqan3::dna5_vector seq{"TTAGGCTGNAG"_dna5};
    int32_t reference_start_position{1};

    auto alignment =
        seqan3::alignment_from_cigar(this->cigar_with_hard_clipping, this->ref, reference_start_position, seq);

    EXPECT_RANGE_EQ(std::get<0>(alignment), this->cigar_with_hard_clipping_gapped_ref);
    EXPECT_RANGE_EQ(std::get<1>(alignment), this->cigar_with_hard_clipping_gapped_seq);
}

TEST_F(alignment_from_cigar, short_cigar_string_with_softclipping)
{
    seqan3::dna5_vector seq = "AGAGGGGGATAACCA"_dna5;

    { // soft clipping at front
        std::vector<seqan3::cigar> short_cigar{{5, 'S'_cigar_operation}, {10, 'M'_cigar_operation}}; // 5S 10M

        auto alignment = seqan3::alignment_from_cigar(short_cigar, this->ref, 0, seq);

        EXPECT_RANGE_EQ(std::get<1>(alignment), "GGGATAACCA"_dna5);
    }

    { // soft clipping at back
        std::vector<seqan3::cigar> short_cigar{{10, 'M'_cigar_operation}, {5, 'S'_cigar_operation}}; // 10M 5S

        auto alignment = seqan3::alignment_from_cigar(short_cigar, this->ref, 0, seq);

        EXPECT_RANGE_EQ(std::get<1>(alignment), "AGAGGGGGAT"_dna5);
    }
}
