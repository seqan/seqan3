// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

TEST(align_config_max_error, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::max_error>));
}

TEST(align_config_max_error, configuration)
{
    {
        seqan3::align_cfg::max_error elem{10};
        seqan3::configuration cfg{elem};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(std::get<seqan3::align_cfg::max_error>(cfg).value)>,
                                  uint32_t>), true);

        EXPECT_EQ(std::get<seqan3::align_cfg::max_error>(cfg).value, 10u);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::max_error{10}};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(std::get<seqan3::align_cfg::max_error>(cfg).value)>,
                                  uint32_t>), true);

        EXPECT_EQ(std::get<seqan3::align_cfg::max_error>(cfg).value, 10u);
    }
}
