// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "fixture/alignment_fixture.hpp"

using namespace seqan3::test::alignment::fixture;

template <auto _fixture>
struct pairwise_alignment_fixture : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture_collection{*_fixture}) const &
    {
        return *_fixture;
    }
};

template <typename fixture_t>
class pairwise_alignment_collection_callback_test : public fixture_t
{};

TYPED_TEST_SUITE_P(pairwise_alignment_collection_callback_test);

TYPED_TEST_P(pairwise_alignment_collection_callback_test, ids)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_sequence2_id{}
        | seqan3::align_cfg::on_result{[&](auto && result)
                                       {
                                           EXPECT_EQ(result.sequence1_id(), result.sequence2_id());
                                       }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, score)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::on_result{[&](auto && result)
                                       {
                                           auto id = result.sequence1_id();
                                           EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                                       }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, end_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::output_end_position{}
        | seqan3::align_cfg::on_result{
            [&](auto && result)
            {
                auto id = result.sequence1_id();
                EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                EXPECT_EQ(result.sequence1_end_position(), fixture.get_end_positions()[id].first);
                EXPECT_EQ(result.sequence2_end_position(), fixture.get_end_positions()[id].second);
            }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, begin_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_begin_position{}
        | seqan3::align_cfg::on_result{
            [&](auto && result)
            {
                auto id = result.sequence1_id();
                EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                EXPECT_EQ(result.sequence1_end_position(), fixture.get_end_positions()[id].first);
                EXPECT_EQ(result.sequence2_end_position(), fixture.get_end_positions()[id].second);
                EXPECT_EQ(result.sequence1_begin_position(), fixture.get_begin_positions()[id].first);
                EXPECT_EQ(result.sequence2_begin_position(), fixture.get_begin_positions()[id].second);
            }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_score{}
        | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_begin_position{}
        | seqan3::align_cfg::output_alignment{}
        | seqan3::align_cfg::on_result{
            [&](auto && result)
            {
                auto id = result.sequence1_id();
                EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                EXPECT_EQ(result.sequence1_end_position(), fixture.get_end_positions()[id].first);
                EXPECT_EQ(result.sequence2_end_position(), fixture.get_end_positions()[id].second);
                EXPECT_EQ(result.sequence1_begin_position(), fixture.get_begin_positions()[id].first);
                EXPECT_EQ(result.sequence2_begin_position(), fixture.get_begin_positions()[id].second);

                auto [gapped_database, gapped_query] = result.alignment();
                EXPECT_RANGE_EQ(gapped_database | seqan3::views::to_char, fixture.get_aligned_sequences1()[id]);
                EXPECT_RANGE_EQ(gapped_query | seqan3::views::to_char, fixture.get_aligned_sequences2()[id]);
            }};

    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_collection_callback_test,
                            ids,
                            score,
                            end_positions,
                            begin_positions,
                            alignment);
