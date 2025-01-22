// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/core/configuration/configuration.hpp>

TEST(band_fixed_size, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::band_fixed_size>));
}

TEST(band_fixed_size, construct)
{
    using seqan3::get;

    constexpr int32_t minus_infinity = std::numeric_limits<int32_t>::lowest();
    constexpr int32_t plus_infinity = std::numeric_limits<int32_t>::max();

    { // Default construct
        seqan3::align_cfg::band_fixed_size band_config{};
        EXPECT_EQ(band_config.lower_diagonal, minus_infinity);
        EXPECT_EQ(band_config.upper_diagonal, plus_infinity);
    }

    { // Construct with parameter
        seqan3::align_cfg::band_fixed_size band_config{seqan3::align_cfg::lower_diagonal{-5},
                                                       seqan3::align_cfg::upper_diagonal{5}};

        EXPECT_EQ(band_config.lower_diagonal, -5);
        EXPECT_EQ(band_config.upper_diagonal, 5);
    }
}

TEST(band_fixed_size, assign)
{
    seqan3::align_cfg::band_fixed_size band_config{};

    band_config.lower_diagonal = -5;
    band_config.upper_diagonal = 5;

    EXPECT_EQ(band_config.lower_diagonal, -5);
    EXPECT_EQ(band_config.upper_diagonal, 5);
}

TEST(band_fixed_size, get_and_assign)
{
    using seqan3::get;

    seqan3::align_cfg::band_fixed_size band_config{seqan3::align_cfg::lower_diagonal{-5},
                                                   seqan3::align_cfg::upper_diagonal{5}};
    seqan3::configuration config{band_config};

    auto & selected_band_config = get<seqan3::align_cfg::band_fixed_size>(config);

    EXPECT_EQ(selected_band_config.lower_diagonal, -5);
    EXPECT_EQ(selected_band_config.upper_diagonal, 5);

    selected_band_config.lower_diagonal = -4;
    selected_band_config.upper_diagonal = 8;

    EXPECT_EQ(get<seqan3::align_cfg::band_fixed_size>(config).lower_diagonal, -4);
    EXPECT_EQ(get<seqan3::align_cfg::band_fixed_size>(config).upper_diagonal, 8);
}
