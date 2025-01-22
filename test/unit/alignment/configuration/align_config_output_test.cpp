// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

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
    seqan3::configuration cfg = seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_end_position{}
                              | seqan3::align_cfg::output_begin_position{} | seqan3::align_cfg::output_alignment{}
                              | seqan3::align_cfg::output_sequence1_id{} | seqan3::align_cfg::output_sequence2_id{};

    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_score>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_end_position>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_begin_position>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_alignment>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_sequence1_id>());
    EXPECT_TRUE(cfg.exists<seqan3::align_cfg::output_sequence2_id>());
}
