// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/parallel.hpp>

int main()
{
    // Enable parallel execution of the search algorithm with 8 threads (and allow 1 error of any type).
    seqan3::configuration cfg1 =
        seqan3::search_cfg::parallel{8} | seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    // Alternative solution: assign to the member variable of the parallel configuration
    seqan3::search_cfg::parallel par_cfg{};
    par_cfg.thread_count = 8;
    seqan3::configuration cfg2 = par_cfg | seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    return 0;
}
