// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    auto seq1 = "TTACGTACGGACTAGCTACAACATTACGGACTAC"_dna4;
    auto seq2 = "GGACGACATGACGTACGACTTTACGTACGACTAGC"_dna4;

    // Configure the output:
    auto output_config = seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{}
                       | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

    // Configure the alignment kernel together with the previous output configuration.
    auto config =
        seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                         seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                         seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}}
        | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                              seqan3::mismatch_score{-2}}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{0}, seqan3::align_cfg::extension_score{-4}}
        | output_config;

    for (auto const & res : seqan3::align_pairwise(std::tie(seq1, seq2), config))
    {
        seqan3::debug_stream << "Score: " << res.score() << '\n';
        seqan3::debug_stream << "Begin: (" << res.sequence1_begin_position() << "," << res.sequence2_begin_position()
                             << ")\n";
        seqan3::debug_stream << "End: (" << res.sequence1_end_position() << "," << res.sequence2_end_position()
                             << ")\n";
        seqan3::debug_stream << "Alignment: \n" << res.alignment() << '\n';
    }
}
