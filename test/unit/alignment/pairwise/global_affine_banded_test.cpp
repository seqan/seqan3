// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/global_affine_banded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_global_affine_banded_testing_types = ::testing::Types<
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::banded::dna4_01>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::banded::dna4_same_sequence_upper_diagonal_0>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::global::affine::banded::dna4_same_sequence_lower_diagonal_0>,
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::banded::dna4_small_band>,
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::banded::dna4_single_diagonal>,
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::global::affine::banded::dna4_large_band>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_global_affine_banded,
                               pairwise_alignment_test,
                               pairwise_global_affine_banded_testing_types, );

struct pairwise_global_affine_banded : public ::testing::Test
{
    using fixture_t = decltype(seqan3::test::alignment::fixture::global::affine::banded::dna4_01);

    fixture_t fixture{seqan3::test::alignment::fixture::global::affine::banded::dna4_01};

    auto & band()
    {
        using seqan3::get;
        return get<seqan3::align_cfg::band_fixed_size>(fixture.config);
    }
};

TEST_F(pairwise_global_affine_banded, invalid_band_lower_diagonal_greater_0)
{
    band().lower_diagonal = 1;
    EXPECT_THROW((align_pairwise(std::tie(fixture.sequence1, fixture.sequence2),
                                 fixture.config | seqan3::align_cfg::output_score{})),
                 seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_global_affine_banded, invalid_band_upper_diagonal_smaller_0)
{
    band().lower_diagonal = -4;
    band().upper_diagonal = -1;
    EXPECT_THROW((align_pairwise(std::tie(fixture.sequence1, fixture.sequence2),
                                 fixture.config | seqan3::align_cfg::output_score{})),
                 seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_global_affine_banded, invalid_band_upper_diagonal_smaller_lower_diagonal)
{
    band().upper_diagonal = -6;
    EXPECT_THROW((align_pairwise(std::tie(fixture.sequence1, fixture.sequence2),
                                 fixture.config | seqan3::align_cfg::output_score{})),
                 seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_global_affine_banded, invalid_band_last_cell_not_covered)
{
    band().upper_diagonal = 5;
    auto result_range = align_pairwise(std::tie(fixture.sequence1, fixture.sequence2),
                                       fixture.config | seqan3::align_cfg::output_score{});
    EXPECT_THROW(result_range.begin(), seqan3::invalid_alignment_configuration);
}

TEST(banded_alignment_issue3266_test, wrong_begin_and_end_position)
{
    using namespace seqan3;
    using namespace seqan3::literals;

    auto const configGeneral = align_cfg::scoring_scheme{nucleotide_scoring_scheme{match_score{1}, mismatch_score{-1}}}
                             | align_cfg::gap_cost_affine{align_cfg::open_score{0}, align_cfg::extension_score{-1}}
                             | align_cfg::method_global{align_cfg::free_end_gaps_sequence1_leading{true},
                                                        align_cfg::free_end_gaps_sequence2_leading{true},
                                                        align_cfg::free_end_gaps_sequence1_trailing{true},
                                                        align_cfg::free_end_gaps_sequence2_trailing{true}};

    auto const configBanded =
        configGeneral | align_cfg::band_fixed_size{align_cfg::lower_diagonal{-40}, align_cfg::upper_diagonal{-20}};

    //0         1         2         3         4
    //01234567890123456789012345678901234567890
    std::pair p{"CGTCTA"_dna4, "AAACCCGGGTTTAAACCCGGGTTTCGTGTACCCCCCCCCCC"_dna4};
    //                      CGTCTA

    auto general_results = align_pairwise(p, configGeneral);
    auto general_res = *std::ranges::begin(general_results);
    auto banded_results = align_pairwise(p, configBanded);
    auto banded_res = *std::ranges::begin(banded_results);

    EXPECT_EQ(general_res.score(), banded_res.score());
    EXPECT_EQ(general_res.sequence1_begin_position(), banded_res.sequence1_begin_position());
    EXPECT_EQ(general_res.sequence2_begin_position(), banded_res.sequence2_begin_position());
    EXPECT_EQ(general_res.sequence1_end_position(), banded_res.sequence1_end_position());
    EXPECT_EQ(general_res.sequence2_end_position(), banded_res.sequence2_end_position());
}
