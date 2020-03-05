// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_edit.hpp>

TEST(align_cfg_edit, is_global)
{
    using type_of_cfg = decltype(seqan3::get<seqan3::align_cfg::mode>(seqan3::align_cfg::edit).value);
    EXPECT_TRUE((std::is_same_v<std::remove_reference_t<type_of_cfg>, seqan3::detail::global_alignment_type>));
}

TEST(align_cfg_edit, is_hamming)
{
    auto scheme = seqan3::get<seqan3::align_cfg::scoring>(seqan3::align_cfg::edit).value;

    for (unsigned i = 0; i < decltype(scheme)::matrix_size; ++i)
    {
        for (unsigned j = 0; j < decltype(scheme)::matrix_size; ++j)
        {
            if (i == j)
                EXPECT_EQ((scheme.score(seqan3::assign_rank_to(i, seqan3::dna15{}),
                                        seqan3::assign_rank_to(j, seqan3::dna15{}))), 0);
            else
                EXPECT_EQ((scheme.score(seqan3::assign_rank_to(i, seqan3::dna15{}),
                                        seqan3::assign_rank_to(j, seqan3::dna15{}))), -1);
        }
    }
}

TEST(align_cfg_edit, is_simple_gap)
{
    auto scheme = seqan3::get<seqan3::align_cfg::gap>(seqan3::align_cfg::edit).value;
    EXPECT_EQ(scheme.get_gap_score(), -1);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
}
