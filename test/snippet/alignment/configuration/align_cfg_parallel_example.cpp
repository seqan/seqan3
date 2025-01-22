// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <thread>

#include <seqan3/alignment/configuration/align_config_parallel.hpp>

int main()
{
    // Enables parallel computation with two threads.
    seqan3::align_cfg::parallel cfg_2{2};

    // Enables parallel computation with the number of concurrent threads supported by the current architecture.
    seqan3::align_cfg::parallel cfg_n{std::thread::hardware_concurrency()};
    ;
}
