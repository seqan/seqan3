// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_result_type.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "fixture/alignment_fixture.hpp"

using namespace seqan3::test::alignment::fixture;

template <
    bool compute_score_matrix_,
    typename database_t,
    typename query_t,
    typename align_cfg_t,
    typename word_t,
    typename is_semi_global_t,
    typename traits_t =
        seqan3::detail::default_edit_distance_trait_type<database_t, query_t, align_cfg_t, is_semi_global_t, word_t>>
struct edit_traits_type : traits_t
{
    static constexpr bool compute_score_matrix = compute_score_matrix_;
    static constexpr bool compute_matrix = compute_score_matrix || traits_t::compute_trace_matrix;
};

template <auto _fixture, typename word_t>
struct global_fixture : public ::testing::Test
{
    auto fixture() -> decltype(alignment_fixture{*_fixture}) const &
    {
        return *_fixture;
    }
    template <bool compute_score_matrix, typename database_t, typename query_t, typename align_cfg_t>
    using edit_traits_type =
        ::edit_traits_type<compute_score_matrix, database_t, query_t, align_cfg_t, word_t, std::false_type>;
};

template <auto _fixture, typename word_t>
struct semi_global_fixture : public global_fixture<_fixture, word_t>
{
    template <bool compute_score_matrix, typename database_t, typename query_t, typename align_cfg_t>
    using edit_traits_type =
        ::edit_traits_type<compute_score_matrix, database_t, query_t, align_cfg_t, word_t, std::true_type>;
};

template <typename fixture_t>
class edit_distance_unbanded_test : public fixture_t
{};

TYPED_TEST_SUITE_P(edit_distance_unbanded_test);

template <template <bool, typename...> typename edit_traits_type,
          bool compute_score_matrix = false,
          typename database_t,
          typename query_t,
          typename align_cfg_t>
auto edit_distance(database_t && database, query_t && query, align_cfg_t && align_cfg)
{
    using align_cfg_with_score_type_t = decltype(align_cfg);

    using alignment_result_value_t = typename seqan3::detail::
        align_result_selector<database_t, query_t, std::remove_reference_t<align_cfg_with_score_type_t>>::type;
    using alignment_result_t = seqan3::alignment_result<alignment_result_value_t>;
    auto align_cfg_with_result_type = align_cfg | seqan3::align_cfg::detail::result_type<alignment_result_t>{};
    using align_config_with_result_type_t = decltype(align_cfg_with_result_type);
    using edit_traits = edit_traits_type<compute_score_matrix, database_t, query_t, align_config_with_result_type_t>;
    using algorithm_t =
        seqan3::detail::edit_distance_unbanded<database_t, query_t, align_config_with_result_type_t, edit_traits>;
    auto alignment = algorithm_t{database, query, align_cfg_with_result_type, edit_traits{}};

    // compute alignment
    alignment(0u, [](auto &&) { /* do nothing */ });
    return alignment;
}

TYPED_TEST_P(edit_distance_unbanded_test, score)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam::template edit_traits_type>(database, query, align_cfg);
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded_test, score_matrix)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_alignment{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam::template edit_traits_type, true>(database, query, align_cfg);
    auto score_matrix = alignment.score_matrix();

    EXPECT_EQ(score_matrix.cols(), database.size() + 1);
    EXPECT_EQ(score_matrix.rows(), query.size() + 1);
    EXPECT_EQ(score_matrix, fixture.score_matrix());
    EXPECT_EQ(alignment.score(), fixture.score);
}

TYPED_TEST_P(edit_distance_unbanded_test, trace_matrix)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_alignment{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam::template edit_traits_type>(database, query, align_cfg);
    auto trace_matrix = alignment.trace_matrix();

    EXPECT_EQ(trace_matrix.cols(), database.size() + 1);
    EXPECT_EQ(trace_matrix.rows(), query.size() + 1);
    EXPECT_EQ(trace_matrix, fixture.trace_matrix());
}

TYPED_TEST_P(edit_distance_unbanded_test, end_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_end_position{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam::template edit_traits_type>(database, query, align_cfg);
    auto end_positions = alignment.end_positions();

    EXPECT_EQ(end_positions.first, fixture.sequence1_end_position);
    EXPECT_EQ(end_positions.second, fixture.sequence2_end_position);
}

TYPED_TEST_P(edit_distance_unbanded_test, begin_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_begin_position{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam::template edit_traits_type>(database, query, align_cfg);
    auto begin_positions = alignment.begin_positions();

    EXPECT_EQ(begin_positions.first, fixture.sequence1_begin_position);
    EXPECT_EQ(begin_positions.second, fixture.sequence2_begin_position);
}

TYPED_TEST_P(edit_distance_unbanded_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_alignment{};

    std::vector database = fixture.sequence1;
    std::vector query = fixture.sequence2;

    auto alignment = edit_distance<TypeParam::template edit_traits_type>(database, query, align_cfg);

    auto && [gapped_database, gapped_query] = alignment.alignment();
    EXPECT_RANGE_EQ(gapped_database | seqan3::views::to_char, fixture.aligned_sequence1);
    EXPECT_RANGE_EQ(gapped_query | seqan3::views::to_char, fixture.aligned_sequence2);
}

REGISTER_TYPED_TEST_SUITE_P(edit_distance_unbanded_test,
                            score,
                            score_matrix,
                            trace_matrix,
                            end_positions,
                            begin_positions,
                            alignment);
