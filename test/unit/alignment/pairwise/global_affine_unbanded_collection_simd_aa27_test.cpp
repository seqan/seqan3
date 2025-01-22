// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "pairwise_alignment_collection_test_template.hpp"

namespace seqan3::test::alignment::collection::simd::global::affine::unbanded
{

static auto aa27_all_same = []()
{
    auto base_fixture = fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10;
    using fixture_t = decltype(base_fixture);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 100; ++i)
        data.push_back(base_fixture);

    return alignment_fixture_collection{base_fixture.config | seqan3::align_cfg::vectorised{}, data};
}();

static auto aa27_different_lengths = []()
{
    auto base_fixture_01 = fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10;
    auto base_fixture_02 = fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_small;
    auto base_fixture_03 = fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_empty_first;
    auto base_fixture_04 = fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_empty_second;
    auto base_fixture_05 = fixture::global::affine::unbanded::aa27_blosum62_gap_1_open_10_empty_both;

    using fixture_t = decltype(base_fixture_01);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 25; ++i)
    {
        data.push_back(base_fixture_01);
        data.push_back(base_fixture_02);
        data.push_back(base_fixture_03);
        data.push_back(base_fixture_04);
        data.push_back(base_fixture_05);
    }

    return alignment_fixture_collection{base_fixture_01.config | seqan3::align_cfg::vectorised{}, data};
}();

} // namespace seqan3::test::alignment::collection::simd::global::affine::unbanded

using pairwise_collection_simd_global_affine_unbanded_testing_types = ::testing::Types<
    pairwise_alignment_fixture<&seqan3::test::alignment::collection::simd::global::affine::unbanded::aa27_all_same>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::collection::simd::global::affine::unbanded::aa27_different_lengths>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_collection_simd_global_affine_unbanded_aa27,
                               pairwise_alignment_collection_test,
                               pairwise_collection_simd_global_affine_unbanded_testing_types, );
