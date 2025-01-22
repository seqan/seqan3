// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_output.hpp>

int main()
{
    // Compute only the end position of the aligned sequences.
    seqan3::configuration cfg = seqan3::align_cfg::output_end_position{};
}
