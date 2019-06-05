// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/range/view/to_char.hpp>

#include "fixture/alignment_fixture.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::test::alignment::fixture;

template <auto _fixture>
struct pairwise_alignment_fixture : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }
};

template <typename fixture_t>
class pairwise_alignment_test : public fixture_t
{};

TYPED_TEST_CASE_P(pairwise_alignment_test);

TYPED_TEST_P(pairwise_alignment_test, score)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_score};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
}

TYPED_TEST_P(pairwise_alignment_test, back_coordinate)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_back_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
}

TYPED_TEST_P(pairwise_alignment_test, front_coordinate)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_front_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
    EXPECT_EQ(res.front_coordinate(), fixture.front_coordinate);
}

TYPED_TEST_P(pairwise_alignment_test, alignment)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_alignment};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.back_coordinate(), fixture.back_coordinate);
    EXPECT_EQ(res.front_coordinate(), fixture.front_coordinate);

    auto && [gapped_database, gapped_query] = res.alignment();
    EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.aligned_sequence1);
    EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.aligned_sequence2);
}

REGISTER_TYPED_TEST_CASE_P(pairwise_alignment_test,
                           score,
                           back_coordinate,
                           front_coordinate,
                           alignment);
