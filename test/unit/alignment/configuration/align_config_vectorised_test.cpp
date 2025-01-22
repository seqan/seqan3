// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_vectorised.hpp>
#include <seqan3/core/configuration/configuration.hpp>

TEST(align_config_vectorised, config_element)
{
    seqan3::configuration cfg{seqan3::align_cfg::vectorised{}};
    EXPECT_TRUE(decltype(cfg)::template exists<seqan3::align_cfg::vectorised>());
}
