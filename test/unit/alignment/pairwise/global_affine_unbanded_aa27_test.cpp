// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_global_affine_unbanded_testing_types = ::testing::Types<
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_small>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_empty_first>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_empty_second>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_empty_both>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_global_affine_unbanded_aa27,
                               pairwise_alignment_test,
                               pairwise_global_affine_unbanded_testing_types, );
