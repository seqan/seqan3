// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <functional>
#include <optional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

// ---------------------------------------------------------------------------------------------------------------------
// individual tests
// ---------------------------------------------------------------------------------------------------------------------

TEST(align_config_parallel, config_element_specialisation)
{
    EXPECT_TRUE((seqan3::config_element_specialisation<seqan3::align_cfg::parallel>));
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
