// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_min_score.hpp>
#include <seqan3/core/configuration/configuration.hpp>

TEST(align_config_min_score, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::min_score>));
}

TEST(align_config_min_score, configuration)
{
    {
        seqan3::align_cfg::min_score elem{-10};
        seqan3::configuration cfg{elem};
        auto min_score = std::get<seqan3::align_cfg::min_score>(cfg);
        EXPECT_TRUE((std::is_same_v<decltype(min_score.score), int32_t>));

        EXPECT_EQ(std::get<seqan3::align_cfg::min_score>(cfg).score, -10);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::min_score{-10}};
        auto min_score = std::get<seqan3::align_cfg::min_score>(cfg);
        EXPECT_TRUE((std::is_same_v<decltype(min_score.score), int32_t>));

        EXPECT_EQ(std::get<seqan3::align_cfg::min_score>(cfg).score, -10);
    }
}
