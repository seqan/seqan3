// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

template <typename test_t>
struct align_cfg_output_test : public ::testing::Test
{};

using test_types = ::testing::Types<seqan3::align_cfg::output_score,
                                    seqan3::align_cfg::output_end_position,
                                    seqan3::align_cfg::output_begin_position,
                                    seqan3::align_cfg::output_alignment,
                                    seqan3::align_cfg::output_sequence1_id,
                                    seqan3::align_cfg::output_sequence2_id>;

INSTANTIATE_TYPED_TEST_SUITE_P(output, pipeable_config_element_test, test_types, );

TEST(align_config_output, score)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::output_score{})>,
                              seqan3::align_cfg::output_score>));
}

TEST(align_config_output, end_position)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::output_end_position{})>,
                              seqan3::align_cfg::output_end_position>));
}

TEST(align_config_output, begin_position)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::output_begin_position{})>,
                              seqan3::align_cfg::output_begin_position>));
}

TEST(align_config_output, alignment)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::output_alignment{})>,
                              seqan3::align_cfg::output_alignment>));
}

TEST(align_config_output, sequence1_id)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::output_sequence1_id{})>,
                              seqan3::align_cfg::output_sequence1_id>));
}

TEST(align_config_output, sequence2_id)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::output_sequence2_id{})>,
                              seqan3::align_cfg::output_sequence2_id>));
}

TEST(align_config_output, combine_outputs)
{
    seqan3::configuration cfg = seqan3::align_cfg::output_score{} |
                                seqan3::align_cfg::output_end_position{} |
                                seqan3::align_cfg::output_begin_position{} |
                                seqan3::align_cfg::output_alignment{} |
                                seqan3::align_cfg::output_sequence1_id{} |
                                seqan3::align_cfg::output_sequence2_id{};


    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_score>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_end_position>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_begin_position>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_alignment>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_sequence1_id>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_sequence2_id>());
}
