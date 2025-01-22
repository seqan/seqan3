// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>

TEST(align_cfg_edit, is_hamming)
{
    auto scheme = seqan3::get<seqan3::align_cfg::scoring_scheme>(seqan3::align_cfg::edit_scheme).scheme;

    constexpr uint32_t dna15_size_v = seqan3::alphabet_size<seqan3::dna15>;

    for (uint32_t i = 0; i < dna15_size_v; ++i)
    {
        for (uint32_t j = 0; j < dna15_size_v; ++j)
        {
            if (i == j)
                EXPECT_EQ((scheme.score(seqan3::assign_rank_to(i, seqan3::dna15{}),
                                        seqan3::assign_rank_to(j, seqan3::dna15{}))),
                          0);
            else
                EXPECT_EQ((scheme.score(seqan3::assign_rank_to(i, seqan3::dna15{}),
                                        seqan3::assign_rank_to(j, seqan3::dna15{}))),
                          -1);
        }
    }
}

TEST(align_cfg_edit, is_simple_gap)
{
    using seqan3::get;

    auto gap_cost_scheme = get<seqan3::align_cfg::gap_cost_affine>(seqan3::align_cfg::edit_scheme);
    EXPECT_EQ(gap_cost_scheme.extension_score, -1);
    EXPECT_EQ(gap_cost_scheme.open_score, 0);
}
