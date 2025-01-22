// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
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

    std::vector vec{"ACGTGACTGACT"_dna4, "ACGAAGACCGAT"_dna4, "ACGTGACTGACT"_dna4, "AGGTACGAGCGACACT"_dna4};

    // Configure the alignment kernel.
    auto config = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::min_score{-7}
                | seqan3::align_cfg::output_score{};

    auto alignment_results = seqan3::align_pairwise(seqan3::views::pairwise_combine(vec), config);
    auto filter_v = std::views::filter(
        [](auto && res)
        {
            return res.score() >= -6;
        });

    for (auto const & result : alignment_results | filter_v)
    {
        seqan3::debug_stream << "Score: " << result.score() << '\n';
    }
}
