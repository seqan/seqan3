// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <concepts>
#include <functional>

#include <seqan3/alignment/configuration/align_config_on_result.hpp>

// -----------------------------------------------------------------------------
// Test capturing various callbacks
// -----------------------------------------------------------------------------

TEST(align_cfg_on_result, with_captureless_lambda)
{
    seqan3::align_cfg::on_result on_result_cfg{[](auto && result)
                                               {
                                                   return result;
                                               }};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ((std::invoke(on_result_cfg.callback, 10)), 10);
}
TEST(align_cfg_on_result, with_capturing_lambda)
{
    int global_result = 0;
    seqan3::align_cfg::on_result on_result_cfg{[&](auto && result)
                                               {
                                                   global_result = result;
                                               }};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ(global_result, 0);
    std::invoke(on_result_cfg.callback, 10);
    EXPECT_EQ(global_result, 10);
}

int my_free_function(int v)
{
    return v;
}

TEST(align_cfg_on_result, with_free_function)
{
    seqan3::align_cfg::on_result on_result_cfg{my_free_function};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ((std::invoke(on_result_cfg.callback, 10)), 10);
}

struct my_function_object
{
    template <typename t>
    t operator()(t && v)
    {
        return std::forward<t>(v);
    }
};

TEST(align_cfg_on_result, with_function_object)
{
    seqan3::align_cfg::on_result on_result_cfg{my_function_object{}};

    EXPECT_TRUE((std::invocable<decltype(on_result_cfg.callback), int>));
    EXPECT_EQ((std::invoke(on_result_cfg.callback, 10)), 10);
}
