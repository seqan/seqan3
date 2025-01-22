// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/configuration/configuration.hpp>

TEST(method_global, access_member_variables)
{
    seqan3::align_cfg::method_global opt{}; // default construction

    // all default to false
    EXPECT_FALSE(opt.free_end_gaps_sequence1_leading);
    EXPECT_FALSE(opt.free_end_gaps_sequence2_leading);
    EXPECT_FALSE(opt.free_end_gaps_sequence1_trailing);
    EXPECT_FALSE(opt.free_end_gaps_sequence2_trailing);

    opt.free_end_gaps_sequence1_leading = true;
    opt.free_end_gaps_sequence2_leading = true;
    opt.free_end_gaps_sequence1_trailing = true;
    opt.free_end_gaps_sequence2_trailing = true;

    EXPECT_TRUE(opt.free_end_gaps_sequence1_leading);
    EXPECT_TRUE(opt.free_end_gaps_sequence2_leading);
    EXPECT_TRUE(opt.free_end_gaps_sequence1_trailing);
    EXPECT_TRUE(opt.free_end_gaps_sequence2_trailing);
}
