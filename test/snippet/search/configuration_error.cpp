// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/max_error.hpp>

int main()
{
    // Allow 1 error of any type.
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};

    // Do not allow substitutions. Allow at most 1 error.
    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                     | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                     | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                     | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}};

    // Sets total errors to 2.
    seqan3::configuration const cfg3 = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                     | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                     | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}};

    // Allow 10% errors of any type.
    seqan3::configuration const cfg4 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{0.1}};

    // Do not allow substitutions. Allow at most 10% errors.
    seqan3::configuration const cfg5 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{.1}}
                                     | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{.0}}
                                     | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{.1}}
                                     | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{.1}};

    // Sets total errors to 20%.
    seqan3::configuration const cfg6 = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{.0}}
                                     | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{.1}}
                                     | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{.1}};

    // Mixed error rate & count: Allow 2 insertions and or 2 deletions and 20% errors in total.

    seqan3::configuration const cfg7 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{.2}}
                                     | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                     | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{2}}
                                     | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{2}};

    return 0;
}
