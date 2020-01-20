// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_edit.hpp>

using namespace seqan3;

TEST(align_cfg_edit, is_global)
{
    EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::mode>(align_cfg::edit).value)>,
                              detail::global_alignment_type>),
               true);
}

TEST(align_cfg_edit, is_hamming)
{
    auto scheme = get<align_cfg::scoring>(align_cfg::edit).value;

    for (unsigned i = 0; i < decltype(scheme)::matrix_size; ++i)
    {
        for (unsigned j = 0; j < decltype(scheme)::matrix_size; ++j)
        {
            if (i == j)
                EXPECT_EQ((scheme.score(assign_rank_to(i, dna15{}), assign_rank_to(j, dna15{}))), 0);
            else
                EXPECT_EQ((scheme.score(assign_rank_to(i, dna15{}), assign_rank_to(j, dna15{}))), -1);
        }
    }
}

TEST(align_cfg_edit, is_simple_gap)
{
    auto scheme = get<align_cfg::gap>(align_cfg::edit).value;
    EXPECT_EQ(scheme.get_gap_score(), -1);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
}
