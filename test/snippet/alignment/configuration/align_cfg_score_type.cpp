// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/core/configuration/configuration.hpp>

int main()
{
    // Compute only the score.
    seqan3::configuration cfg1 =
        seqan3::align_cfg::score_type<int16_t>{};                        // Now the alignment computes 16 bit integers.
    seqan3::configuration cfg2 = seqan3::align_cfg::score_type<float>{}; // Now the alignment computes float scores.
}
