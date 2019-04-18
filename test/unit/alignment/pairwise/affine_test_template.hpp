// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include <seqan3/range/view/to_char.hpp>

#include "alignment_fixture.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::test::alignment::fixture;

template <auto _fixture>
struct param : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }
};

template <typename param_t>
class align_affine : public param_t
{};

TYPED_TEST_CASE_P(align_affine);

TYPED_TEST_P(align_affine, score)
{
    auto const & fixture = this->fixture();
    // We only compute the score.
    auto align_cfg = fixture.config | align_cfg::result{with_score};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    // TODO Make test work with ranges.
    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    EXPECT_EQ((*std::ranges::begin(alignment)).score(), fixture.score);
}

TYPED_TEST_P(align_affine, end_position)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_back_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    auto res = *std::ranges::begin(alignment);
    EXPECT_EQ(res.score(), fixture.score);

    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
}

TYPED_TEST_P(align_affine, begin_position)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_front_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    auto res = *std::ranges::begin(alignment);
    EXPECT_EQ(res.score(), fixture.score);

    EXPECT_EQ(res.front_coordinate(), fixture.front_coordinate);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
}

TYPED_TEST_P(align_affine, trace)
{
    auto const & fixture = this->fixture();
    auto align_cfg = fixture.config | align_cfg::result{with_alignment};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = align_pairwise(std::tie(database, query), align_cfg);

    auto res = *std::ranges::begin(alignment);
    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.front_coordinate(), fixture.front_coordinate);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
    EXPECT_TRUE(ranges::equal(get<0>(res.alignment()) | view::to_char, fixture.aligned_sequence1));
    EXPECT_TRUE(ranges::equal(get<1>(res.alignment()) | view::to_char, fixture.aligned_sequence2));
}

REGISTER_TYPED_TEST_CASE_P(align_affine, score, end_position, begin_position, trace);
