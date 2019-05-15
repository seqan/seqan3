// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>

namespace seqan3::test::alignment::fixture
{

static constexpr auto INF = std::nullopt;

static constexpr auto NON = detail::trace_directions::none;
static constexpr auto D = detail::trace_directions::diagonal;
static constexpr auto U = detail::trace_directions::up;
static constexpr auto L = detail::trace_directions::left;
static constexpr auto DU = D | U;
static constexpr auto UL = U | L;
static constexpr auto DL = D | L;
static constexpr auto DUL = D | U | L;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t,
          typename score_vector_or_matrix_t, typename trace_vector_or_matrix_t>
struct alignment_fixture
{
    sequence1_t sequence1;
    sequence2_t sequence2;

    config_t config;

    score_t score;
    std::string aligned_sequence1;
    std::string aligned_sequence2;

    alignment_coordinate front_coordinate;
    alignment_coordinate back_coordinate;

    score_vector_or_matrix_t score_vector{};
    trace_vector_or_matrix_t trace_vector{};

    auto score_matrix() const requires seqan3::detail::Matrix<score_vector_or_matrix_t>
    {
        return detail::debug_matrix{score_vector};
    };

    auto score_matrix() const requires !seqan3::detail::Matrix<score_vector_or_matrix_t>
    {
        detail::row_wise_matrix score_matrix{score_vector, sequence2.size() + 1, sequence1.size() + 1};
        return detail::debug_matrix{std::move(score_matrix)};
    };

    auto trace_matrix() const requires seqan3::detail::Matrix<trace_vector_or_matrix_t>
    {
        return detail::debug_matrix{trace_vector};
    };

    auto trace_matrix() const requires !seqan3::detail::Matrix<trace_vector_or_matrix_t>
    {
        detail::row_wise_matrix trace_matrix{trace_vector, sequence2.size() + 1, sequence1.size() + 1};
        return detail::debug_matrix{std::move(trace_matrix)};
    };
};

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t,
          typename score_vector_or_matrix_t, typename trace_vector_or_matrix_t>
alignment_fixture(
    sequence1_t sequence1,
    sequence2_t sequence2,
    config_t config,
    score_t score,
    std::string aligned_sequence1,
    std::string aligned_sequence2,
    alignment_coordinate front_coordinate,
    alignment_coordinate back_coordinate,
    score_vector_or_matrix_t score_vector,
    trace_vector_or_matrix_t trace_vector
)
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t, score_vector_or_matrix_t, trace_vector_or_matrix_t>;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t>
alignment_fixture(
    sequence1_t sequence1,
    sequence2_t sequence2,
    config_t config,
    score_t score,
    std::string aligned_sequence1,
    std::string aligned_sequence2,
    alignment_coordinate front_coordinate,
    alignment_coordinate back_coordinate
)
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t,
                     std::vector<score_t>, std::vector<detail::trace_directions>>;


} // namespace seqan3::test::alignment::fixture
