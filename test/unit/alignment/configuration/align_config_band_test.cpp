// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

TEST(align_config_band, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::band<seqan3::static_band>>));
}

TEST(align_config_band, configuration)
{
    {
        seqan3::align_cfg::band elem{seqan3::static_band{seqan3::lower_bound{-5}, seqan3::upper_bound{5}}};
        seqan3::configuration cfg{elem};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::band>(cfg).value)>,
                                  seqan3::static_band>), true);

        EXPECT_EQ(seqan3::get<seqan3::align_cfg::band>(cfg).value.lower_bound, -5);
        EXPECT_EQ(seqan3::get<seqan3::align_cfg::band>(cfg).value.upper_bound, 5);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-5},
                                                                              seqan3::upper_bound{5}}}};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::band>(cfg).value)>,
                                  seqan3::static_band>), true);

        EXPECT_EQ(seqan3::get<seqan3::align_cfg::band>(cfg).value.lower_bound, -5);
        EXPECT_EQ(seqan3::get<seqan3::align_cfg::band>(cfg).value.upper_bound, 5);
    }
}
