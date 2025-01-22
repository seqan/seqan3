// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/output.hpp>

int main()
{
    auto const zero_errors = seqan3::search_cfg::error_count{0};
    // No errors, all hits as text position
    seqan3::configuration const default_cfg =
        seqan3::search_cfg::max_error_total{zero_errors} | seqan3::search_cfg::max_error_substitution{zero_errors}
        | seqan3::search_cfg::max_error_insertion{zero_errors} | seqan3::search_cfg::max_error_deletion{zero_errors}
        | seqan3::search_cfg::output_query_id{} | seqan3::search_cfg::output_reference_id{}
        | seqan3::search_cfg::output_reference_begin_position{} | seqan3::search_cfg::hit_all{};
    return 0;
}
