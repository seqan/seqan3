// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/configuration/configuration.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// pipeable_config_element_test template
// ---------------------------------------------------------------------------------------------------------------------

using config_element_types = ::testing::Types<seqan3::align_cfg::method_global,
                                              seqan3::align_cfg::method_local>;

INSTANTIATE_TYPED_TEST_SUITE_P(method, pipeable_config_element_test, config_element_types, );

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
