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
#include <seqan3/utility/range/to.hpp>
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
class pairwise_alignment_collection_test : public fixture_t
{};

TYPED_TEST_SUITE_P(pairwise_alignment_collection_test);

TYPED_TEST_P(pairwise_alignment_collection_test, score)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_score{};
    auto [database, query] = fixture.get_sequences();
    auto alignment_rng = seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);

    EXPECT_RANGE_EQ(alignment_rng
                        | std::views::transform(
                            [](auto res)
                            {
                                return res.score();
                            }),
                    fixture.get_scores());
}

TYPED_TEST_P(pairwise_alignment_collection_test, end_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg =
        fixture.config | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_end_position{};

    using traits_t = seqan3::detail::alignment_configuration_traits<decltype(align_cfg)>;

    if constexpr (!(traits_t::is_vectorised && traits_t::is_banded))
    {
        auto [database, query] = fixture.get_sequences();
        auto res_vec =
            seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg) | seqan3::ranges::to<std::vector>();

        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return res.score();
                                }),
                        fixture.get_scores());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::pair<size_t, size_t>{res.sequence1_end_position(),
                                                                     res.sequence2_end_position()};
                                }),
                        fixture.get_end_positions());
    }
}

TYPED_TEST_P(pairwise_alignment_collection_test, begin_positions)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_begin_position{}
                                    | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{};

    using traits_t = seqan3::detail::alignment_configuration_traits<decltype(align_cfg)>;

    if constexpr (!traits_t::is_vectorised)
    {
        auto [database, query] = fixture.get_sequences();
        auto res_vec =
            seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg) | seqan3::ranges::to<std::vector>();

        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return res.score();
                                }),
                        fixture.get_scores());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::pair<size_t, size_t>{res.sequence1_end_position(),
                                                                     res.sequence2_end_position()};
                                }),
                        fixture.get_end_positions());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::pair<size_t, size_t>{res.sequence1_begin_position(),
                                                                     res.sequence2_begin_position()};
                                }),
                        fixture.get_begin_positions());
    }
}

TYPED_TEST_P(pairwise_alignment_collection_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::output_alignment{}
                                    | seqan3::align_cfg::output_begin_position{}
                                    | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_score{};

    using traits_t = seqan3::detail::alignment_configuration_traits<decltype(align_cfg)>;

    if constexpr (!traits_t::is_vectorised)
    {
        auto [database, query] = fixture.get_sequences();
        auto res_vec =
            seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg) | seqan3::ranges::to<std::vector>();

        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return res.score();
                                }),
                        fixture.get_scores());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::pair<size_t, size_t>{res.sequence1_end_position(),
                                                                     res.sequence2_end_position()};
                                }),
                        fixture.get_end_positions());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::pair<size_t, size_t>{res.sequence1_begin_position(),
                                                                     res.sequence2_begin_position()};
                                }),
                        fixture.get_begin_positions());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::get<0>(res.alignment()) | seqan3::views::to_char
                                         | seqan3::ranges::to<std::string>();
                                }),
                        fixture.get_aligned_sequences1());
        EXPECT_RANGE_EQ(res_vec
                            | std::views::transform(
                                [](auto res)
                                {
                                    return std::get<1>(res.alignment()) | seqan3::views::to_char
                                         | seqan3::ranges::to<std::string>();
                                }),
                        fixture.get_aligned_sequences2());
    }
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_collection_test, score, end_positions, begin_positions, alignment);
