// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/std/concepts>

TEST(align_config_gap, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::gap<seqan3::gap_scheme<>>>));
}

TEST(align_config_gap, configuration)
{
    {
        seqan3::align_cfg::gap elem{seqan3::gap_scheme<>{}};
        seqan3::configuration cfg{elem};
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::gap>(cfg).value.get_gap_score()), -1);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::gap{seqan3::gap_scheme<>{}}};
        EXPECT_EQ((seqan3::get<seqan3::align_cfg::gap>(cfg).value.get_gap_open_score()), 0);
    }
}
