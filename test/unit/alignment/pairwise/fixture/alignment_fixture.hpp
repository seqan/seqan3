// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>

namespace seqan3::test::alignment::fixture
{

static constexpr auto INF = detail::matrix_inf<int>;

static constexpr auto NON = detail::trace_directions::none;
static constexpr auto D = detail::trace_directions::diagonal;
static constexpr auto U = detail::trace_directions::up;
static constexpr auto L = detail::trace_directions::left;
static constexpr auto DU = D | U;
static constexpr auto UL = U | L;
static constexpr auto DL = D | L;
static constexpr auto DUL = D | U | L;

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename trace_t>
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

    std::vector<score_t> score_vector{};
    std::vector<trace_t> trace_vector{};

    detail::alignment_score_matrix<std::vector<score_t>> score_matrix() const
    {
        return {score_vector, sequence2.size() + 1, sequence1.size() + 1};
    };

    detail::alignment_trace_matrix<std::vector<trace_t>> trace_matrix() const
    {
        return {trace_vector, sequence2.size() + 1, sequence1.size() + 1};
    };
};

template <typename sequence1_t, typename sequence2_t, typename config_t, typename score_t, typename trace_t>
alignment_fixture(
    sequence1_t sequence1,
    sequence2_t sequence2,
    config_t config,
    score_t score,
    std::string aligned_sequence1,
    std::string aligned_sequence2,
    alignment_coordinate front_coordinate,
    alignment_coordinate back_coordinate,
    std::vector<score_t> score_vector,
    std::vector<trace_t> trace_vector
)
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t, trace_t>;

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
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t, detail::trace_directions>;

template <typename config_t, typename alignment_fixture_t>
struct alignment_fixture_collection
{
    config_t config;
    std::vector<alignment_fixture_t> collection;

    auto get_sequences() const
    {
        std::vector<decltype(collection[0].sequence1)> vec1;
        std::vector<decltype(collection[0].sequence2)> vec2;

        for (size_t i = 0; i < collection.size(); ++i)
        {
            vec1.push_back(collection[i].sequence1);
            vec2.push_back(collection[i].sequence2);
        }
        return std::pair{vec1, vec2};
    }

    auto get_scores() const
    {
        std::vector<decltype(collection[0].score)> vec;
        for (size_t i = 0; i < collection.size(); ++i)
            vec.push_back(collection[i].score);

        return vec;
    }

    auto get_back_coordinates() const
    {
        std::vector<decltype(collection[0].back_coordinate)> vec;
        for (size_t i = 0; i < collection.size(); ++i)
            vec.push_back(collection[i].back_coordinate);

        return vec;
    }

    auto get_front_coordinates() const
    {
        std::vector<decltype(collection[0].front_coordinate)> vec;
        for (size_t i = 0; i < collection.size(); ++i)
            vec.push_back(collection[i].front_coordinate);

        return vec;
    }

    auto get_aligned_sequences1() const
    {
        std::vector<decltype(collection[0].aligned_sequence1)> vec;
        for (size_t i = 0; i < collection.size(); ++i)
            vec.push_back(collection[i].aligned_sequence1);

        return vec;
    }

    auto get_aligned_sequences2() const
    {
        std::vector<decltype(collection[0].aligned_sequence2)> vec;
        for (size_t i = 0; i < collection.size(); ++i)
            vec.push_back(collection[i].aligned_sequence2);

        return vec;
    }
};

template <typename config_t, typename fixture_t>
alignment_fixture_collection(
    config_t config,
    std::vector<fixture_t> collection
)
-> alignment_fixture_collection<config_t, fixture_t>;

} // namespace seqan3::test::alignment::fixture
