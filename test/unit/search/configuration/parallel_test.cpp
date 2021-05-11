// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <optional>
#include <seqan3/search/configuration/parallel.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// test template : pipeable_config_element_test
// ---------------------------------------------------------------------------------------------------------------------

using test_types = ::testing::Types<seqan3::search_cfg::parallel>;

INSTANTIATE_TYPED_TEST_SUITE_P(parallel_elements, pipeable_config_element_test, test_types, );

// ---------------------------------------------------------------------------------------------------------------------
// individual tests
// ---------------------------------------------------------------------------------------------------------------------

TEST(search_config_parallel, member_variable)
{
    {   // default construction
        seqan3::search_cfg::parallel cfg{};
        EXPECT_FALSE(cfg.thread_count);
        EXPECT_THROW(cfg.thread_count.value(), std::bad_optional_access);
    }

    {   // construct with value
        seqan3::search_cfg::parallel cfg{4};
        EXPECT_EQ(cfg.thread_count.value(), 4u);
    }

    {   // assign value
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
