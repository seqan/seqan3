// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>


#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include <seqan3/range/view/to_char.hpp>

#include "fixture/global_affine_unbanded.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::fixture;

template <auto _fixture>
struct param : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }
};

template <typename param_t>
class global_affine_unbanded : public param_t
{};

TYPED_TEST_CASE_P(global_affine_unbanded);

using global_affine_unbanded_types
    = ::testing::Types<
        param<&global::affine::unbanded::dna4_01>
    >;

TYPED_TEST_P(global_affine_unbanded, score)
{
    auto const & fixture = this->fixture();
    // We only compute the score.
    auto align_cfg = fixture.config | align_cfg::result{align_cfg::with_score};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    // TODO Make test work with ranges.
    auto alignment = align_pairwise(std::pair{database, query}, align_cfg);

    EXPECT_EQ((*std::ranges::begin(alignment)).get_score(), fixture.score);
}

REGISTER_TYPED_TEST_CASE_P(global_affine_unbanded, score);

// work around a bug that you can't specify more than 50 template arguments to ::testing::types
INSTANTIATE_TYPED_TEST_CASE_P(global, global_affine_unbanded, global_affine_unbanded_types);
