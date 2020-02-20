// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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

static constexpr auto N     = seqan3::detail::trace_directions::none;
static constexpr auto D     = seqan3::detail::trace_directions::diagonal;
static constexpr auto u     = seqan3::detail::trace_directions::up;
static constexpr auto l     = seqan3::detail::trace_directions::left;
static constexpr auto U     = seqan3::detail::trace_directions::up_open;
static constexpr auto L     = seqan3::detail::trace_directions::left_open;
static constexpr auto DU    = D | U;
static constexpr auto UL    = U | L;
static constexpr auto DL    = D | L;
static constexpr auto DUL   = D | U | L;
static constexpr auto Du    = D | u;
static constexpr auto Uu    = U | u;
static constexpr auto uL    = u | L;
static constexpr auto DuL   = D | u | L;
static constexpr auto DUu   = D | U | u;
static constexpr auto UuL   = U | u | L;
static constexpr auto DUuL  = D | U | u | L;
static constexpr auto Dl    = D | l;
static constexpr auto Ul    = U | l;
static constexpr auto DUl   = D | U | l;
static constexpr auto Ll    = L | l;
static constexpr auto DLl   = D | L | l;
static constexpr auto ULl   = U | L | l;
static constexpr auto DULl  = D | U | L | l;
static constexpr auto ul    = u | l;
static constexpr auto Dul   = D | u | l;
static constexpr auto Uul   = U | u | l ;
static constexpr auto DUul  = D | U | u | l;
static constexpr auto uLl   = u | L | l;
static constexpr auto DuLl  = D | u | L | l;
static constexpr auto UuLl  = U | u | L | l;
static constexpr auto DUuLl = D | U | u | L | l;

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

    seqan3::alignment_coordinate front_coordinate;
    seqan3::alignment_coordinate back_coordinate;

    score_vector_or_matrix_t score_vector{};
    trace_vector_or_matrix_t trace_vector{};

    auto score_matrix() const requires seqan3::detail::matrix<score_vector_or_matrix_t>
    {
        return seqan3::detail::debug_matrix{score_vector};
    };

    auto score_matrix() const requires !seqan3::detail::matrix<score_vector_or_matrix_t>
    {
        seqan3::detail::row_wise_matrix<std::ranges::range_value_t<score_vector_or_matrix_t>>
            score_matrix{seqan3::detail::number_rows{sequence2.size() + 1},
                         seqan3::detail::number_cols{sequence1.size() + 1},
                         score_vector};
        return seqan3::detail::debug_matrix{std::move(score_matrix)};
    };

    auto trace_matrix() const requires seqan3::detail::matrix<trace_vector_or_matrix_t>
    {
        return seqan3::detail::debug_matrix{trace_vector};
    };

    auto trace_matrix() const requires !seqan3::detail::matrix<trace_vector_or_matrix_t>
    {
        seqan3::detail::row_wise_matrix<std::ranges::range_value_t<trace_vector_or_matrix_t>>
            trace_matrix{seqan3::detail::number_rows{sequence2.size() + 1},
                         seqan3::detail::number_cols{sequence1.size() + 1},
                         trace_vector};
        return seqan3::detail::debug_matrix{std::move(trace_matrix)};
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
    seqan3::alignment_coordinate front_coordinate,
    seqan3::alignment_coordinate back_coordinate,
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
    seqan3::alignment_coordinate front_coordinate,
    seqan3::alignment_coordinate back_coordinate
)
-> alignment_fixture<sequence1_t, sequence2_t, config_t, score_t,
                     std::vector<score_t>, std::vector<seqan3::detail::trace_directions>>;

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
