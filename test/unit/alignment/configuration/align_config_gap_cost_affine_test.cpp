// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/concepts>

#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/core/configuration/configuration.hpp>

TEST(align_config_gap, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::gap_cost_affine>));
}

TEST(align_config_gap, configuration)
{
    using seqan3::get;
    {
        seqan3::align_cfg::gap_cost_affine elem{}; // default construction
        seqan3::configuration cfg{elem};
        EXPECT_EQ((get<seqan3::align_cfg::gap_cost_affine>(cfg).extension_score), -1);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::gap_cost_affine{}}; // default construction
        EXPECT_EQ((get<seqan3::align_cfg::gap_cost_affine>(cfg).open_score), 0);
    }

    {
        seqan3::align_cfg::gap_cost_affine scheme{seqan3::align_cfg::open_score{-10},
                                                  seqan3::align_cfg::extension_score{-1}};
        EXPECT_EQ((scheme.open_score), -10);
        EXPECT_EQ((scheme.extension_score), -1);
    }

    {
        seqan3::align_cfg::gap_cost_affine scheme{seqan3::align_cfg::open_score{-10},
                                                  seqan3::align_cfg::extension_score{-1}};
        EXPECT_EQ((scheme.open_score), -10);
        EXPECT_EQ((scheme.extension_score), -1);
    }
}
