// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "fixture/alignment_fixture.hpp"

template <auto _fixture>
struct pairwise_alignment_fixture : public ::testing::Test
{
    auto fixture() -> decltype(seqan3::test::alignment::fixture::alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }
};

template <typename fixture_t>
class pairwise_alignment_test : public fixture_t
{};

TYPED_TEST_SUITE_P(pairwise_alignment_test);

TYPED_TEST_P(pairwise_alignment_test, score)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
}

TYPED_TEST_P(pairwise_alignment_test, end_positions)
{
    auto const & fixture = this->fixture();

    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::result{seqan3::with_score,
                                                                                 seqan3::using_score_type<double>} |
                                                       seqan3::align_cfg::output_score |
                                                       seqan3::align_cfg::output_end_position;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_TRUE((std::same_as<decltype(res.score()), double>));
    EXPECT_EQ(res.sequence1_end_position(), fixture.end_positions.first);
    EXPECT_EQ(res.sequence2_end_position(), fixture.end_positions.second);
}

TYPED_TEST_P(pairwise_alignment_test, begin_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score |
                                                       seqan3::align_cfg::output_end_position |
                                                       seqan3::align_cfg::output_begin_position;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.sequence1_end_position(), fixture.end_positions.first);
    EXPECT_EQ(res.sequence2_end_position(), fixture.end_positions.second);
    EXPECT_EQ(res.sequence1_begin_position(), fixture.begin_positions.first);
    EXPECT_EQ(res.sequence2_begin_position(), fixture.begin_positions.second);
}

TYPED_TEST_P(pairwise_alignment_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score |
                                                       seqan3::align_cfg::output_end_position |
                                                       seqan3::align_cfg::output_begin_position |
                                                       seqan3::align_cfg::output_alignment |
                                                       seqan3::align_cfg::detail::debug;

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.sequence1_end_position(), fixture.end_positions.first);
    EXPECT_EQ(res.sequence2_end_position(), fixture.end_positions.second);
    EXPECT_EQ(res.sequence1_begin_position(), fixture.begin_positions.first);
    EXPECT_EQ(res.sequence2_begin_position(), fixture.begin_positions.second);

    auto && [gapped_database, gapped_query] = res.alignment();
    EXPECT_EQ(gapped_database | seqan3::views::to_char | seqan3::views::to<std::string>, fixture.aligned_sequence1);
    EXPECT_EQ(gapped_query | seqan3::views::to_char | seqan3::views::to<std::string>, fixture.aligned_sequence2);

    using score_matrix_t = seqan3::detail::two_dimensional_matrix<std::optional<int32_t>>;
    using trace_matrix_t = seqan3::detail::two_dimensional_matrix<std::optional<seqan3::detail::trace_directions>>;

    EXPECT_RANGE_EQ(static_cast<score_matrix_t>(res.score_matrix()), fixture.score_vector);
    EXPECT_RANGE_EQ(static_cast<trace_matrix_t>(res.trace_matrix()), fixture.trace_vector);
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_test, score, end_positions, begin_positions, alignment);
