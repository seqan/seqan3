// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_global_affine_unbanded_testing_types = ::testing::Types<
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_01>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_02>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_03>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_04>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_05>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_seq1_empty>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_seq2_empty>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_both_empty>,
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::issue_3043>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_global_affine_unbanded,
                               pairwise_alignment_test,
                               pairwise_global_affine_unbanded_testing_types, );
