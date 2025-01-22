// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::get;

    seqan3::configuration my_cfg =
        seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
        | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4},
                                             seqan3::align_cfg::upper_diagonal{4}};
    // my_cfg is now of type configuration<gap_cost_affine, band_fixed_size>

    seqan3::debug_stream << get<1>(my_cfg).lower_diagonal << '\n';                                   // prints -4
    seqan3::debug_stream << get<seqan3::align_cfg::band_fixed_size>(my_cfg).upper_diagonal << '\n';  // prints 4
    seqan3::debug_stream << get<seqan3::align_cfg::gap_cost_affine>(my_cfg).extension_score << '\n'; // prints -1
}
