// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <functional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_vectorise.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

TEST(align_config_vectorise, config_element)
{
    seqan3::configuration cfg{seqan3::align_cfg::vectorise};
    EXPECT_TRUE(decltype(cfg)::template exists<seqan3::detail::vectorise_tag>());
}
