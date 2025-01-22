// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    // Basic alignment algorithm configuration.
    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

    std::pair p{"ACGTAGC"_dna4, "AGTACGACG"_dna4};

    // Compute only the score:
    for (auto res : seqan3::align_pairwise(p, config | seqan3::align_cfg::output_score{}))
        seqan3::debug_stream << res << "\n"; // prints: {score: -4}

    // Compute only the alignment:
    for (auto res : seqan3::align_pairwise(p, config | seqan3::align_cfg::output_alignment{}))
        seqan3::debug_stream << res << "\n"; // prints: {alignment: (ACGTA-G-C-,A-GTACGACG)}

    // Compute the score and the alignment:
    for (auto res :
         seqan3::align_pairwise(p, config | seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_alignment{}))
        seqan3::debug_stream << res << "\n"; // prints: {score: -4, alignment: (ACGTA-G-C-,A-GTACGACG)}

    // By default compute everything:
    for (auto res : seqan3::align_pairwise(p, config))
        seqan3::debug_stream
            << res << "\n"; // prints {id: 0, score: -4, begin: (0,0), end: (7,9) alignment: (ACGTA-G-C-,A-GTACGACG)}
}
