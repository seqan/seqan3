// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//![start]
#include <vector>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;
    //![start]

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

    {
        //![example1]
        std::pair p{"ACGTAGC"_dna4, "AGTACGACG"_dna4};
        auto result = seqan3::align_pairwise(p, config);
        //![example1]
    }

    {
        //![example2]
        std::vector vec{"ACCA"_dna4, "ATTA"_dna4};
        auto result = seqan3::align_pairwise(std::tie(vec[0], vec[1]), config);
        //![example2]
    }

    //![example3]
    std::vector vec{std::pair{"AGTGCTACG"_dna4, "ACGTGCGACTAG"_dna4},
                    std::pair{"AGTAGACTACG"_dna4, "ACGTACGACACG"_dna4},
                    std::pair{"AGTTACGAC"_dna4, "AGTAGCGATCG"_dna4}};

    // Compute the alignment of a single pair.
    auto edit_config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;

    for (auto const & res : seqan3::align_pairwise(std::tie(vec[0].first, vec[0].second), edit_config))
        seqan3::debug_stream << "The score: " << res.score() << "\n";

    // Compute the alignment over a range of pairs.
    for (auto const & res : seqan3::align_pairwise(vec, edit_config))
        seqan3::debug_stream << "The score: " << res.score() << "\n";
    //![example3]
    //![end]
}
//![end]
