// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_output.hpp>

int main()
{
    // Compute only the score.
    seqan3::configuration cfg = seqan3::align_cfg::output_score{};
}
