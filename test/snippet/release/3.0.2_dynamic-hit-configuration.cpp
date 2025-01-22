// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

int main()
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                          "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                          "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};
    seqan3::fm_index index{text};

    // Use the dynamic hit configuration to set hit_all_best mode.
    seqan3::configuration search_config = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit{seqan3::search_cfg::hit_all_best{}};

    seqan3::debug_stream << "All single best hits:\n";
    for (auto && hit : search(query, index, search_config)) // Find all best hits:
        seqan3::debug_stream << hit << '\n';

    // Change the hit configuration to the strata mode with a stratum of 1.
    using seqan3::get;
    get<seqan3::search_cfg::hit>(search_config).hit_variant = seqan3::search_cfg::hit_strata{1};

    seqan3::debug_stream << "\nAll x+1 hits:\n";
    for (auto && hit : search(query, index, search_config)) // Find all strata hits.
        seqan3::debug_stream << hit << '\n';
}
