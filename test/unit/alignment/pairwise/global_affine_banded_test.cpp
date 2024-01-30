// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
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
