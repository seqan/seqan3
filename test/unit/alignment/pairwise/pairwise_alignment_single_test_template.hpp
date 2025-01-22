// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>

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
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
}

TYPED_TEST_P(pairwise_alignment_test, end_positions)
{
    auto const & fixture = this->fixture();

    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_end_position{}
                                    | seqan3::align_cfg::output_score{} | seqan3::align_cfg::score_type<double>{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_SAME_TYPE(decltype(res.score()), double);
    EXPECT_EQ(res.sequence1_end_position(), fixture.sequence1_end_position);
    EXPECT_EQ(res.sequence2_end_position(), fixture.sequence2_end_position);
}

TYPED_TEST_P(pairwise_alignment_test, begin_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_begin_position{}
                                    | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.sequence1_end_position(), fixture.sequence1_end_position);
    EXPECT_EQ(res.sequence2_end_position(), fixture.sequence2_end_position);
    EXPECT_EQ(res.sequence1_begin_position(), fixture.sequence1_begin_position);
    EXPECT_EQ(res.sequence2_begin_position(), fixture.sequence2_begin_position);
}

TYPED_TEST_P(pairwise_alignment_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score{}
                                    | seqan3::align_cfg::output_end_position{}
                                    | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_alignment{}
                                    | seqan3::align_cfg::detail::debug{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment_rng = seqan3::align_pairwise(std::tie(database, query), align_cfg);
    auto res = *alignment_rng.begin();

    EXPECT_EQ(res.score(), fixture.score);
    EXPECT_EQ(res.sequence1_end_position(), fixture.sequence1_end_position);
    EXPECT_EQ(res.sequence2_end_position(), fixture.sequence2_end_position);
    EXPECT_EQ(res.sequence1_begin_position(), fixture.sequence1_begin_position);
    EXPECT_EQ(res.sequence2_begin_position(), fixture.sequence2_begin_position);

    auto && [gapped_database, gapped_query] = res.alignment();
    EXPECT_RANGE_EQ(gapped_database | seqan3::views::to_char, fixture.aligned_sequence1);
    EXPECT_RANGE_EQ(gapped_query | seqan3::views::to_char, fixture.aligned_sequence2);

    using score_matrix_t = seqan3::detail::two_dimensional_matrix<std::optional<int32_t>>;
    using trace_matrix_t = seqan3::detail::two_dimensional_matrix<std::optional<seqan3::detail::trace_directions>>;

    EXPECT_RANGE_EQ(static_cast<score_matrix_t>(res.score_matrix()), fixture.score_vector);
    EXPECT_RANGE_EQ(static_cast<trace_matrix_t>(res.trace_matrix()), fixture.trace_vector);
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_test, score, end_positions, begin_positions, alignment);
