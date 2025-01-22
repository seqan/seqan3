// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3::literals;

int main()
{
    // configure a global alignment for DNA sequences
    auto min_cfg = seqan3::align_cfg::method_global{}
                 | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                       seqan3::mismatch_score{-5}}}
                 | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                      seqan3::align_cfg::extension_score{-1}};

    auto seq1 = "TCGT"_dna4;
    auto seq2 = "ACGA"_dna4;
    for (auto res : seqan3::align_pairwise(std::tie(seq1, seq2), min_cfg))
        seqan3::debug_stream << res.score() << '\n'; // print out the alignment score
}
