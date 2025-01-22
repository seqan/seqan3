// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_min_score.hpp>

int main()
{
    // Computes semi global edit distance using fast-bit vector algorithm.
    auto cfg_fast = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
                  | seqan3::align_cfg::edit_scheme;

    // Computes semi global edit distance using slower standard pairwise algorithm.
    auto cfg_slow = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
                  | seqan3::align_cfg::edit_scheme;

    // Computes global distance allowing a minimal score of 3 (Default: edit distance).
    auto cfg_errors =
        seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-3};
}
