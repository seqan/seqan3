// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/cigar_conversion/cigar_from_alignment.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
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

struct cigar_from_alignment : public cigar_conversion_data
{};

TEST_F(cigar_from_alignment, empty_sequences)
{
    std::vector<seqan3::gapped<seqan3::dna5>> empty{};
    EXPECT_THROW(seqan3::cigar_from_alignment(std::tie(empty, empty)), std::logic_error);
}

TEST_F(cigar_from_alignment, aligned_sequences_do_not_have_the_same_length)
{
    std::vector<seqan3::gapped<seqan3::dna5>> too_short{'A'_dna5};
    EXPECT_THROW(seqan3::cigar_from_alignment(std::tie(this->simple_cigar_gapped_ref, too_short)), std::logic_error);
}

TEST_F(cigar_from_alignment, simple_cigar)
{
    auto cigar = seqan3::cigar_from_alignment(std::tie(this->simple_cigar_gapped_ref, this->simple_cigar_gapped_seq),
                                              {.soft_front = 1});

    EXPECT_RANGE_EQ(cigar, this->simple_cigar);
}

TEST_F(cigar_from_alignment, with_padding)
{
    auto cigar =
        seqan3::cigar_from_alignment(std::tie(this->cigar_with_padding_gapped_ref, this->cigar_with_padding_gapped_seq),
                                     {.soft_front = 1, .soft_back = 1});

    EXPECT_RANGE_EQ(cigar, this->cigar_with_padding);
}

TEST_F(cigar_from_alignment, extended_cigar)
{
    auto cigar =
        seqan3::cigar_from_alignment(std::tie(this->cigar_with_padding_gapped_ref, this->cigar_with_padding_gapped_seq),
                                     {.soft_front = 1, .soft_back = 1},
                                     true /* output extended cigar */);

    EXPECT_RANGE_EQ(cigar, this->extended_cigar_with_padding);
}

TEST_F(cigar_from_alignment, hard_clipping)
{
    auto cigar = seqan3::cigar_from_alignment(
        std::tie(this->cigar_with_hard_clipping_gapped_ref, this->cigar_with_hard_clipping_gapped_seq),
        {.hard_front = 1, .hard_back = 2, .soft_front = 2, .soft_back = 1});

    EXPECT_RANGE_EQ(cigar, this->cigar_with_hard_clipping);
}
