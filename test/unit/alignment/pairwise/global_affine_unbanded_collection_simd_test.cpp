// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "pairwise_alignment_collection_test_template.hpp"

namespace seqan3::test::alignment::collection::simd::global::affine::unbanded
{

static auto dna4_all_same = []()
{
    auto base_fixture = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_01;
    using fixture_t = decltype(base_fixture);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 100; ++i)
        data.push_back(base_fixture);

    return alignment_fixture_collection{base_fixture.config | seqan3::align_cfg::vectorise, data};
}();

static auto dna4_different_length = []()
{
    auto base_fixture_01 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_01;
    auto base_fixture_02 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_02;
    auto base_fixture_03 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_03;
    auto base_fixture_04 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_04;

    using fixture_t = decltype(base_fixture_01);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 25; ++i)
    {
        data.push_back(base_fixture_01);
        data.push_back(base_fixture_02);
        data.push_back(base_fixture_03);
        data.push_back(base_fixture_04);
    }

    return alignment_fixture_collection{base_fixture_01.config | seqan3::align_cfg::vectorise, data};
}();

static auto dna4_with_empty_sequences = []()
{
    auto base_fixture_01 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_01;
    auto base_fixture_02 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_02;
    auto base_fixture_03 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_part_04;
    auto base_fixture_04 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_seq1_empty;
    auto base_fixture_05 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_seq2_empty;
    auto base_fixture_06 = fixture::global::affine::unbanded::dna4_match_4_mismatch_5_gap_1_open_10_both_empty;

    using fixture_t = decltype(base_fixture_01);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 25; ++i)
    {
        data.push_back(base_fixture_01);
        data.push_back(base_fixture_02);
        data.push_back(base_fixture_03);
        data.push_back(base_fixture_04);
        data.push_back(base_fixture_05);
        data.push_back(base_fixture_06);
    }

    return alignment_fixture_collection{base_fixture_01.config | seqan3::align_cfg::vectorise, data};
}();

} // namespace seqan3::test::alignment::collection::simd::global::affine::unbanded

using pairwise_collection_simd_global_affine_unbanded_testing_types = ::testing::Types<
        pairwise_alignment_fixture<&seqan3::test::alignment::collection::simd::global::affine::unbanded::dna4_all_same>,
        pairwise_alignment_fixture<&seqan3::test::alignment::collection::simd::global::affine::unbanded::dna4_different_length>,
        pairwise_alignment_fixture<&seqan3::test::alignment::collection::simd::global::affine::unbanded::dna4_with_empty_sequences>
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_collection_simd_global_affine_unbanded,
                               pairwise_alignment_collection_test,
                               pairwise_collection_simd_global_affine_unbanded_testing_types, );
