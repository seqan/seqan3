// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector text{
        "CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};

    seqan3::fm_index index{text};

    seqan3::debug_stream << "Searching all hits\n";
    seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                              | seqan3::search_cfg::hit{seqan3::search_cfg::hit_all{}};
    auto results_all = search(query, index, cfg);
    // Attention: results_all is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There are " << std::ranges::distance(results_all) << " hits.\n";

    seqan3::debug_stream << "Searching all best hits\n";
    using seqan3::get;
    get<seqan3::search_cfg::hit>(cfg).hit_variant = seqan3::search_cfg::hit_all_best{};
    auto results_all_best = search(query, index, cfg);
    // Attention: results_all_best is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There are " << std::ranges::distance(results_all_best) << " hits.\n";

    seqan3::debug_stream << "Searching best hit\n";
    get<seqan3::search_cfg::hit>(cfg).hit_variant = seqan3::search_cfg::hit_single_best{};
    auto results_best = search(query, index, cfg);
    // Attention: results_best is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There is " << std::ranges::distance(results_best) << " hits.\n";

    seqan3::debug_stream << "Searching all hits in the 1-stratum\n";
    get<seqan3::search_cfg::hit>(cfg).hit_variant = seqan3::search_cfg::hit_strata{1};
    auto results_strata = search(query, index, cfg);
    // Attention: results_strata is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There are " << std::ranges::distance(results_strata) << " hits.\n";
}
