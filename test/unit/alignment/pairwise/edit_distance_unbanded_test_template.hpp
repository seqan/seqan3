// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>

#include <seqan3/range/view/to_char.hpp>

#include "fixture/alignment_fixture.hpp"

using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::test::alignment::fixture;

template <typename word_t = uint64_t, typename is_semi_global_t = std::false_type>
struct test_traits_type
{
    using word_type = word_t;
    using is_semi_global_type = is_semi_global_t;
};

template <auto _fixture, typename word_t>
struct global_fixture : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }

    using traits_type = test_traits_type<word_t, std::false_type>;
};

template <auto _fixture, typename word_t>
struct semi_global_fixture : public global_fixture<_fixture, word_t>
{
    using traits_type = test_traits_type<word_t, std::true_type>;
};

template <typename fixture_t>
class edit_distance_unbanded : public fixture_t
{};

TYPED_TEST_CASE_P(edit_distance_unbanded);

template <typename traits_t, typename database_t, typename query_t, typename align_cfg_t>
auto edit_distance(database_t && database, query_t && query, align_cfg_t && align_cfg)
{
    using algorithm_t = pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_cfg_t, traits_t>;
    auto alignment = algorithm_t{database, query, align_cfg};

    // compute alignment
    alignment(0u);
    return alignment;
}

TYPED_TEST_P(edit_distance_unbanded, score)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_score};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<typename TypeParam::traits_type>(database, query, align_cfg);
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded, score_matrix)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_alignment};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<typename TypeParam::traits_type>(database, query, align_cfg);
    auto score_matrix = alignment.score_matrix();

    EXPECT_EQ(score_matrix.cols(), database.size()+1);
    EXPECT_EQ(score_matrix.rows(), query.size()+1);
    EXPECT_EQ(score_matrix, fixture.score_matrix());
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded, trace_matrix)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_alignment};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<typename TypeParam::traits_type>(database, query, align_cfg);
    auto trace_matrix = alignment.trace_matrix();

    EXPECT_EQ(trace_matrix.cols(), database.size()+1);
    EXPECT_EQ(trace_matrix.rows(), query.size()+1);
    EXPECT_EQ(trace_matrix, fixture.trace_matrix());
}

TYPED_TEST_P(edit_distance_unbanded, back_coordinate)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_back_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<typename TypeParam::traits_type>(database, query, align_cfg);
    auto back_coordinate = alignment.back_coordinate();

    EXPECT_EQ(back_coordinate, fixture.back_coordinate);
}

TYPED_TEST_P(edit_distance_unbanded, front_coordinate)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_front_coordinate};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<typename TypeParam::traits_type>(database, query, align_cfg);
    auto front_coordinate = alignment.front_coordinate();

    EXPECT_EQ(front_coordinate, fixture.front_coordinate);
}

TYPED_TEST_P(edit_distance_unbanded, alignment)
{
    auto const & fixture = this->fixture();
    configuration align_cfg = fixture.config | align_cfg::result{with_alignment};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<typename TypeParam::traits_type>(database, query, align_cfg);

    auto && [gapped_database, gapped_query] = alignment.alignment();
    EXPECT_EQ(std::string{gapped_database | view::to_char}, fixture.aligned_sequence1);
    EXPECT_EQ(std::string{gapped_query | view::to_char}, fixture.aligned_sequence2);
}

REGISTER_TYPED_TEST_CASE_P(edit_distance_unbanded,
                           score,
                           score_matrix,
                           trace_matrix,
                           back_coordinate,
                           front_coordinate,
                           alignment);
