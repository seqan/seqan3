// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/semi_global_affine_unbanded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_semiglobal_affine_unbanded_testing_types = ::testing::Types<
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::unbanded::dna4_01_semi_first>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::unbanded::dna4_02_semi_first>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::unbanded::dna4_03_semi_second>,
        pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::unbanded::dna4_04_semi_second>
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_semiglobal_affine_unbanded,
                               pairwise_alignment_test,
                               pairwise_semiglobal_affine_unbanded_testing_types, );
