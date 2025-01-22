// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>

int main()
{
    // Report all hits with 0 errors (maximum number of errors defaults to 0).
    seqan3::configuration const cfg1 = seqan3::search_cfg::hit_all{};

    // Report all hits with 0 and 1 errors.
    seqan3::configuration const cfg2 =
        seqan3::search_cfg::hit_all{} | seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    // Report the single best hit with the least number of errors (up to 1 error is allowed).
    seqan3::configuration const cfg3 =
        seqan3::search_cfg::hit_single_best{} | seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    // Report all hits with the least number of errors (either 0 or 1 errors).
    seqan3::configuration const cfg4 =
        seqan3::search_cfg::hit_all_best{} | seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    // Report all hits with best + 1 error but no more than 2 (errors).
    // E.g., if the best hit has 1 error, all hits with 1 and 2 errors are reported.
    // E.g., if the best hit has 2 error, only hits with 2 errors are reported since 3 exceeds total.
    seqan3::configuration const cfg5 =
        seqan3::search_cfg::hit_strata{1} | seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{2}};

    // you must choose only one mode
    // auto fail = seqan3::search_cfg::hit_single_best{} | seqan3::search_cfg::hit_all{}; // doesn't compile

    return 0;
}
