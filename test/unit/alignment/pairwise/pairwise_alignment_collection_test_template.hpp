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
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

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
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::result{seqan3::with_score};
    auto [database, query] = fixture.get_sequences();
    auto alignment_rng = seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg);

    auto scores = fixture.get_scores();
    EXPECT_TRUE((std::ranges::equal(alignment_rng | std::views::transform([] (auto res) { return res.score(); } ),
                                    fixture.get_scores())));
}

TYPED_TEST_P(pairwise_alignment_collection_test, back_coordinate)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::result{seqan3::with_back_coordinate};
    auto [database, query] = fixture.get_sequences();
    auto res_vec = seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg)
                 | seqan3::views::to<std::vector>;

    EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.score(); }),
                                    fixture.get_scores())));
    EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.back_coordinate(); }),
                                    fixture.get_back_coordinates())));
}

TYPED_TEST_P(pairwise_alignment_collection_test, front_coordinate)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::result{seqan3::with_front_coordinate};

    using traits_t = seqan3::detail::alignment_configuration_traits<decltype(align_cfg)>;

    if constexpr (!traits_t::is_vectorised)
    {
        auto [database, query] = fixture.get_sequences();
        auto res_vec = seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg)
                     | seqan3::views::to<std::vector>;

        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.score(); }),
                                        fixture.get_scores())));
        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.back_coordinate(); }),
                                        fixture.get_back_coordinates())));
        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.front_coordinate(); }),
                                        fixture.get_front_coordinates())));
    }
}

TYPED_TEST_P(pairwise_alignment_collection_test, alignment)
{
    auto const & fixture = this->fixture();
    seqan3::configuration align_cfg = fixture.config | seqan3::align_cfg::result{seqan3::with_alignment};

    using traits_t = seqan3::detail::alignment_configuration_traits<decltype(align_cfg)>;

    if constexpr (!traits_t::is_vectorised)
    {
        auto [database, query] = fixture.get_sequences();
        auto res_vec = seqan3::align_pairwise(seqan3::views::zip(database, query), align_cfg)
                     | seqan3::views::to<std::vector>;

        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.score(); }),
                                        fixture.get_scores())));
        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.back_coordinate(); }),
                                        fixture.get_back_coordinates())));
        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res) { return res.front_coordinate(); }),
                                        fixture.get_front_coordinates())));
        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res)
                                                {
                                                    return std::get<0>(res.alignment()) | seqan3::views::to_char
                                                                                        | seqan3::views::to<std::string>;
                                                }),
                                        fixture.get_aligned_sequences1())));
        EXPECT_TRUE((std::ranges::equal(res_vec | std::views::transform([] (auto res)
                                                {
                                                    return std::get<1>(res.alignment()) | seqan3::views::to_char
                                                                                        | seqan3::views::to<std::string>;
                                                }),
                                        fixture.get_aligned_sequences2())));
    }
}

REGISTER_TYPED_TEST_SUITE_P(pairwise_alignment_collection_test, score, back_coordinate, front_coordinate, alignment);
