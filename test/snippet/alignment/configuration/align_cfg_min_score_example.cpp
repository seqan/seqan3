// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_min_score.hpp>
#include <seqan3/core/configuration/configuration.hpp>

int main()
{
    // Allow a minimal score of -5, i.e. at most 5 edit operations.
    seqan3::configuration config = seqan3::align_cfg::min_score{-5};
    auto min_score = std::get<seqan3::align_cfg::min_score>(config);
    min_score.score = -5;
}
