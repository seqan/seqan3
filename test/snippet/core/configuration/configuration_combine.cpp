// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

int main()
{
    // Here we use the global and banded alignment configurations to show how they can be combined.
    seqan3::configuration my_cfg = seqan3::align_cfg::method_global{}
                                 | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4},
                                                                      seqan3::align_cfg::upper_diagonal{4}};
}
