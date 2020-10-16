// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <vector>

#include <gtest/gtest.h>

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
        pairwise_alignment_fixture<&seqan3::test::alignment::collection::simd::global::affine::banded::dna4_all_same>
    >;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_collection_simd_global_affine_banded,
                               pairwise_alignment_collection_test,
                               pairwise_collection_simd_global_affine_banded_testing_types, );
