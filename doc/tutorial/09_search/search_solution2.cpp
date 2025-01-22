// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

using namespace seqan3::literals;

void run_text_single()
{
    seqan3::dna4_vector text{
        "CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "=====   Running on a single text   =====\n"
                         << "The following hits were found:\n";

    for (auto && result : search("GCT"_dna4, index))
        seqan3::debug_stream << result << '\n';
}

void run_text_collection()
{
    std::vector<seqan3::dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                          "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                          "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "===== Running on a text collection =====\n"
                         << "The following hits were found:\n";

    for (auto && result : search("GCT"_dna4, index))
        seqan3::debug_stream << result << '\n';
}

int main()
{
    run_text_single();
    seqan3::debug_stream << '\n';
    run_text_collection();
}
