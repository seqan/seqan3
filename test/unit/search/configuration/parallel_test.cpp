// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <optional>

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/parallel.hpp>

TEST(search_config_parallel, member_variable)
{
    { // default construction
        seqan3::search_cfg::parallel cfg{};
        EXPECT_FALSE(cfg.thread_count);
        EXPECT_THROW(cfg.thread_count.value(), std::bad_optional_access);
    }

    { // construct with value
        seqan3::search_cfg::parallel cfg{4};
        EXPECT_EQ(cfg.thread_count.value(), 4u);
    }

    { // assign value
        seqan3::search_cfg::parallel cfg{};
        cfg.thread_count = 4;
        EXPECT_EQ(cfg.thread_count.value(), 4u);
    }
}

TEST(search_config_parallel, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::search_cfg::parallel>));
}

TEST(search_config_parallel, configuration)
{
    { // from lvalue.
        seqan3::search_cfg::parallel elem{4};
        seqan3::configuration cfg{elem};
        using ret_type = decltype(std::get<seqan3::search_cfg::parallel>(cfg).thread_count);
        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<ret_type>, std::optional<uint32_t>>));

        EXPECT_EQ(std::get<seqan3::search_cfg::parallel>(cfg).thread_count.value(), 4u);
    }

    { // from rvalue.
        seqan3::configuration cfg{seqan3::search_cfg::parallel{4}};
        using ret_type = decltype(std::get<seqan3::search_cfg::parallel>(cfg).thread_count);
        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<ret_type>, std::optional<uint32_t>>));

        EXPECT_EQ(std::get<seqan3::search_cfg::parallel>(cfg).thread_count.value(), 4u);
    }
}
