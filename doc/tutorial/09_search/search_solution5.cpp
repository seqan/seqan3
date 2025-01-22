// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <span>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

using namespace seqan3::literals;

// Define the pairwise alignment configuration globally.
inline constexpr auto align_config =
    seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                     seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                     seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
    | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_alignment{} | seqan3::align_cfg::output_score{};

void run_text_single()
{
    seqan3::dna4_vector text{
        "CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "Searching all best hits allowing for 1 error in a single text\n";

    seqan3::configuration const search_config =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_all_best{};

    auto search_results = search(query, index, search_config);

    seqan3::debug_stream << "-----------------\n";

    for (auto && hit : search_results)
    {
        size_t start = hit.reference_begin_position() ? hit.reference_begin_position() - 1 : 0;
        std::span text_view{std::data(text) + start, query.size() + 1};

        for (auto && res : align_pairwise(std::tie(text_view, query), align_config))
        {
            auto && [aligned_database, aligned_query] = res.alignment();
            seqan3::debug_stream << "score:    " << res.score() << '\n';
            seqan3::debug_stream << "database: " << aligned_database << '\n';
            seqan3::debug_stream << "query:    " << aligned_query << '\n';
            seqan3::debug_stream << "=============\n";
        }
    }
}

void run_text_collection()
{
    std::vector<seqan3::dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                          "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                          "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "Searching all best hits allowing for 1 error in a text collection\n";

    seqan3::configuration const search_config =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_all_best{};

    seqan3::debug_stream << "-----------------\n";

    for (auto && hit : search(query, index, search_config))
    {
        size_t start = hit.reference_begin_position() ? hit.reference_begin_position() - 1 : 0;
        std::span text_view{std::data(text[hit.reference_id()]) + start, query.size() + 1};

        for (auto && res : align_pairwise(std::tie(text_view, query), align_config))
        {
            auto && [aligned_database, aligned_query] = res.alignment();
            seqan3::debug_stream << "score:    " << res.score() << '\n';
            seqan3::debug_stream << "database: " << aligned_database << '\n';
            seqan3::debug_stream << "query:    " << aligned_query << '\n';
            seqan3::debug_stream << "=============\n";
        }
    }
}

int main()
{
    run_text_single();
    seqan3::debug_stream << '\n';
    run_text_collection();
}
