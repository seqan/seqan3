// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/local_affine_banded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_local_affine_banded_testing_types =
    ::testing::Types<pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::dna4_01>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::dna4_02>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::dna4_03>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::dna4_04>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::dna4_05>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::dna4_06>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::rna5_01>,
                     pairwise_alignment_fixture<&seqan3::test::alignment::fixture::local::affine::banded::aa27_01>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_local_affine_banded,
                               pairwise_alignment_test,
                               pairwise_local_affine_banded_testing_types, );
