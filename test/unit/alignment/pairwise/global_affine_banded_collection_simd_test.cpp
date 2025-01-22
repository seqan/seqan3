// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_banded.hpp"
#include "pairwise_alignment_collection_test_template.hpp"

namespace seqan3::test::alignment::collection::simd::global::affine::banded
{

static auto dna4_all_same = []()
{
    auto base_fixture = fixture::global::affine::banded::dna4_01;
    using fixture_t = decltype(base_fixture);

    std::vector<fixture_t> data{};
    for (size_t i = 0; i < 100; ++i)
        data.push_back(base_fixture);

    return alignment_fixture_collection{base_fixture.config | seqan3::align_cfg::vectorised{}, data};
}();

} // namespace seqan3::test::alignment::collection::simd::global::affine::banded

using pairwise_collection_simd_global_affine_banded_testing_types = ::testing::Types<
    pairwise_alignment_fixture<&seqan3::test::alignment::collection::simd::global::affine::banded::dna4_all_same>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_collection_simd_global_affine_banded,
                               pairwise_alignment_collection_test,
                               pairwise_collection_simd_global_affine_banded_testing_types, );
