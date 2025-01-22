// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>

#include "fixture/semi_global_affine_banded.hpp"
#include "pairwise_alignment_single_test_template.hpp"

using pairwise_semiglobal_affine_banded_testing_types = ::testing::Types<
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_01_semi_first>,
    // pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_02_semi_first>,
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_03_semi_second>,
    pairwise_alignment_fixture<&seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_04_semi_second>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_free_lb_with_band_tl2br_no_matches>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_free_tlbr_with_band_tl2b>,
    pairwise_alignment_fixture<
        &seqan3::test::alignment::fixture::semi_global::affine::banded::dna4_free_tlbr_with_band_tl2r>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pairwise_semiglobal_affine_banded,
                               pairwise_alignment_test,
                               pairwise_semiglobal_affine_banded_testing_types, );

using seqan3::operator""_dna4;

struct pairwise_semiglobal_affine_banded : public ::testing::Test
{
    seqan3::dna4_vector sequence1 = "TTTTTACGTATGTCCCCC"_dna4;
    seqan3::dna4_vector sequence2 = "ACGTAAAACGTAGTAGAT"_dna4;

    auto base_config() const
    {
        return seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                  seqan3::align_cfg::extension_score{-1}}
             | seqan3::align_cfg::scoring_scheme{
                 seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};
    }
};

TEST_F(pairwise_semiglobal_affine_banded, invalid_band_lower_diagonal_greater_0)
{
    auto config = base_config()
                | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{1},
                                                     seqan3::align_cfg::upper_diagonal{8}};

    EXPECT_THROW((align_pairwise(std::tie(sequence1, sequence2), config | seqan3::align_cfg::output_score{})),
                 seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_semiglobal_affine_banded, invalid_band_upper_diagonal_smaller_0)
{
    auto config = base_config()
                | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-8},
                                                     seqan3::align_cfg::upper_diagonal{-1}};

    EXPECT_THROW((align_pairwise(std::tie(sequence1, sequence2), config | seqan3::align_cfg::output_score{})),
                 seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_semiglobal_affine_banded, invalid_band_upper_diagonal_smaller_lower_diagonal)
{
    auto config = base_config()
                | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{3},
                                                     seqan3::align_cfg::upper_diagonal{-3}};

    EXPECT_THROW((align_pairwise(std::tie(sequence1, sequence2), config | seqan3::align_cfg::output_score{})),
                 seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_semiglobal_affine_banded, invalid_band_lower_diagonal_ends_in_last_column)
{
    auto config = base_config()
                | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
                | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{4},
                                                     seqan3::align_cfg::upper_diagonal{8}};
    auto alignment_range = align_pairwise(std::tie(sequence1, sequence2), config | seqan3::align_cfg::output_score{});
    EXPECT_THROW((alignment_range.begin()), seqan3::invalid_alignment_configuration);
}

TEST_F(pairwise_semiglobal_affine_banded, invalid_band_upper_diagonal_ends_in_last_row)
{
    auto config = base_config()
                | seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-8},
                                                     seqan3::align_cfg::upper_diagonal{-4}};

    auto alignment_range = align_pairwise(std::tie(sequence1, sequence2), config | seqan3::align_cfg::output_score{});
    EXPECT_THROW((alignment_range.begin()), seqan3::invalid_alignment_configuration);
}
