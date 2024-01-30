// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>

TEST(align_cfg_edit, is_hamming)
{
    auto scheme = seqan3::get<seqan3::align_cfg::scoring_scheme>(seqan3::align_cfg::edit_scheme).scheme;

    for (unsigned i = 0; i < decltype(scheme)::matrix_size; ++i)
    {
        for (unsigned j = 0; j < decltype(scheme)::matrix_size; ++j)
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
