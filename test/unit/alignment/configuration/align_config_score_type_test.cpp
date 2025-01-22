// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

TEST(align_config_score_type, score_type)
{
    EXPECT_TRUE((std::same_as<std::remove_cvref_t<decltype(seqan3::align_cfg::score_type<int32_t>{})>,
                              seqan3::align_cfg::score_type<int32_t>>));                 // default case
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
