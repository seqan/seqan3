// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <functional>
#include <optional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/core/configuration/configuration.hpp>

// ---------------------------------------------------------------------------------------------------------------------
// individual tests
// ---------------------------------------------------------------------------------------------------------------------

TEST(align_config_parallel, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::parallel>));
}

TEST(align_config_parallel, configuration)
{
    { // from lvalue.
        seqan3::align_cfg::parallel elem{2};
        seqan3::configuration cfg{elem};
        auto cfg_value = std::get<seqan3::align_cfg::parallel>(cfg).thread_count;

        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(cfg_value)>, std::optional<uint32_t>>));
        EXPECT_EQ(cfg_value, 2u);
    }

    { // from rvalue.
        seqan3::configuration cfg{seqan3::align_cfg::parallel{2}};
        auto cfg_value = std::get<seqan3::align_cfg::parallel>(cfg).thread_count;

        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(cfg_value)>, std::optional<uint32_t>>));
        EXPECT_EQ(cfg_value, 2u);
    }
}
