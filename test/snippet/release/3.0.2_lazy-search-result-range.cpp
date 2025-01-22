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

    seqan3::configuration const search_config =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_all_best{};

    // Always provide a unified interface over the found hits independent of the index its text layout.
    seqan3::debug_stream << "Search in text collection:\n";
    seqan3::fm_index index_collection{text};
    for (auto && hit : search(query, index_collection, search_config)) // Over a text collection.
        seqan3::debug_stream << hit << '\n';

    seqan3::debug_stream << "\nSearch in single text:\n";
    seqan3::fm_index index_single{text[0]};
    for (auto && hit : search(query, index_single, search_config)) // Over a single text.
        seqan3::debug_stream << hit << '\n';
}
