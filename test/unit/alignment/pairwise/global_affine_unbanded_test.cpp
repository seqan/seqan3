// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_global_affine_unbanded_testing_types = ::testing::Types<
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_01>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_02>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_03>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_04>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_05>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_seq1_empty>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_seq2_empty>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_both_empty>
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_global_affine_unbanded,
                               pairwise_alignment_test,
                               pairwise_global_affine_unbanded_testing_types, );
