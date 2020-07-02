// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/test/expect_range_eq.hpp>

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

TYPED_TEST_P(pairwise_alignment_collection_callback_test, score)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_score} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          auto id = result.id();
                                          EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                                      }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, back_coordinate)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_back_coordinate} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          auto id = result.id();
                                          EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                                          EXPECT_EQ(result.sequence1_end_position(),
                                                    fixture.get_back_coordinates()[id].first);
                                          EXPECT_EQ(result.sequence2_end_position(),
                                                    fixture.get_back_coordinates()[id].second);
                                      }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, front_coordinate)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_front_coordinate} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          auto id = result.id();
                                          EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                                          EXPECT_EQ(result.sequence1_end_position(),
                                                    fixture.get_back_coordinates()[id].first);
                                          EXPECT_EQ(result.sequence2_end_position(),
                                                    fixture.get_back_coordinates()[id].second);
                                          EXPECT_EQ(result.sequence1_begin_position(),
                                                    fixture.get_front_coordinates()[id].first);
                                          EXPECT_EQ(result.sequence2_begin_position(),
                                                    fixture.get_front_coordinates()[id].second);
                                      }};
    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_collection_callback_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_alignment} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          auto id = result.id();
                                          EXPECT_EQ(result.score(), fixture.get_scores()[id]);
                                          EXPECT_EQ(result.sequence1_end_position(),
                                                    fixture.get_back_coordinates()[id].first);
                                          EXPECT_EQ(result.sequence2_end_position(),
                                                    fixture.get_back_coordinates()[id].second);
                                          EXPECT_EQ(result.sequence1_begin_position(),
                                                    fixture.get_front_coordinates()[id].first);
                                          EXPECT_EQ(result.sequence2_begin_position(),
                                                    fixture.get_front_coordinates()[id].second);

                                          auto [gapped_database, gapped_query] = result.alignment();
                                          EXPECT_RANGE_EQ(gapped_database | seqan3::views::to_char,
                                                          fixture.get_aligned_sequences1()[id]);
                                          EXPECT_RANGE_EQ(gapped_query | seqan3::views::to_char,
                                                          fixture.get_aligned_sequences2()[id]);
                                      }};

    auto [database, query] = fixture.get_sequences();
    seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_collection_callback_test,
                            score,
                            back_coordinate,
                            front_coordinate,
                            alignment);
