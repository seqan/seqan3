// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

template <typename test_t>
struct align_cfg_output_test : public ::testing::Test
{};

using test_types = ::testing::Types<seqan3::align_cfg::score_type<int32_t>>;

INSTANTIATE_TYPED_TEST_SUITE_P(score_type, pipeable_config_element_test, test_types, );

TEST(align_config_score_type, score_type)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::score_type<int32_t>{})>,
                              seqan3::align_cfg::score_type<int32_t>>));                            // default case
    EXPECT_SAME_TYPE(decltype(seqan3::align_cfg::score_type<int32_t>{})::type, int32_t); // default case

    EXPECT_SAME_TYPE(decltype(seqan3::align_cfg::score_type<int16_t>{})::type, int16_t);
    EXPECT_SAME_TYPE(decltype(seqan3::align_cfg::score_type<float>{})::type, float);
    EXPECT_SAME_TYPE(decltype(seqan3::align_cfg::score_type<double>{})::type, double);
}

TEST(align_config_score_type, score_type_exists)
{
    seqan3::configuration cfg = seqan3::align_cfg::score_type<double>{};
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::score_type<double>>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::score_type>());
}
