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
class pairwise_alignment_callback_test : public fixture_t
{};

TYPED_TEST_SUITE_P(pairwise_alignment_callback_test);

TYPED_TEST_P(pairwise_alignment_callback_test, score)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_score} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          EXPECT_EQ(result.score(), fixture.score);
                                      }};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;
    seqan3::align_pairwise(std::tie(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_callback_test, back_coordinate)
{
    auto const & fixture = this->fixture();

    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_back_coordinate,
                                                                seqan3::using_score_type<double>} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          EXPECT_EQ(result.score(), fixture.score);
                                          EXPECT_TRUE((std::same_as<decltype(result.score()), double>));
                                          EXPECT_EQ(result.back_coordinate(), fixture.back_coordinate);
                                      }};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    seqan3::align_pairwise(std::tie(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_callback_test, front_coordinate)
{
    auto const & fixture = this->fixture();

    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_front_coordinate} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          EXPECT_EQ(result.score(), fixture.score);
                                          EXPECT_EQ(result.back_coordinate(), fixture.back_coordinate);
                                          EXPECT_EQ(result.front_coordinate(), fixture.front_coordinate);
                                      }};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    seqan3::align_pairwise(std::tie(database, query), align_cfg);
}

TYPED_TEST_P(pairwise_alignment_callback_test, alignment)
{
    auto const & fixture = this->fixture();

    seqan3::configuration align_cfg = fixture.config |
                                      seqan3::align_cfg::result{seqan3::with_alignment} |
                                      seqan3::align_cfg::on_result{[&] (auto && result)
                                      {
                                          EXPECT_EQ(result.score(), fixture.score);
                                          EXPECT_EQ(result.back_coordinate(), fixture.back_coordinate);
                                          EXPECT_EQ(result.front_coordinate(), fixture.front_coordinate);

                                          auto && [gapped_database, gapped_query] = result.alignment();
                                          EXPECT_RANGE_EQ(gapped_database | seqan3::views::to_char,
                                                          fixture.aligned_sequence1);
                                          EXPECT_RANGE_EQ(gapped_query | seqan3::views::to_char,
                                                          fixture.aligned_sequence2);
                                      }};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    seqan3::align_pairwise(std::tie(database, query), align_cfg);
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_callback_test,
                            score,
                            back_coordinate,
                            front_coordinate,
                            alignment);
