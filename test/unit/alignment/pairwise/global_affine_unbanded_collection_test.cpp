// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_unbanded.hpp"
#include "pairwise_alignment_collection_test_template.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::test::alignment;

namespace seqan3::test::alignment::collection::global::affine::unbanded
{

static auto dna4_01 = []()
{

    using fixture_t = decltype(fixture::global::affine::unbanded::dna4_01);
    std::vector<fixture_t> data;
    for (size_t i = 0; i < 100; ++i)
        data.push_back(fixture::global::affine::unbanded::dna4_01);

    return alignment_fixture_collection{fixture::global::affine::unbanded::dna4_01.config, data};
}();

} // namespace seqan3::test::alignment::fixture::collection::global::affine::unbanded

using pairwise_global_affine_collection_unbanded_testing_types = ::testing::Types<
        pairwise_alignment_fixture<&collection::global::affine::unbanded::dna4_01, sequenced_policy>,
        pairwise_alignment_fixture<&collection::global::affine::unbanded::dna4_01, parallel_policy>
    >;

INSTANTIATE_TYPED_TEST_CASE_P(pairwise_global_affine_collection_unbanded,
                              pairwise_alignment_collection_test,
                              pairwise_global_affine_collection_unbanded_testing_types);
