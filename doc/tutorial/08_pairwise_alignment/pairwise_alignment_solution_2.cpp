// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector vec{"ACGTGAACTGACT"_dna4, "ACGAAGACCGAT"_dna4, "ACGTGACTGACT"_dna4, "AGGTACGAGCGACACT"_dna4};

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
                | seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{}};

    for (auto const & res : seqan3::align_pairwise(seqan3::views::pairwise_combine(vec), config))
        seqan3::debug_stream << "Score: " << res.score() << '\n';
}
