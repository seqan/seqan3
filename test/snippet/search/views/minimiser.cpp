// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"ACGTAGC"_dna4};

    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{3}});
    seqan3::debug_stream << hashes << '\n'; // [6,27,44,50,9]

    auto minimiser = hashes | seqan3::views::minimiser(4);
    seqan3::debug_stream << minimiser << '\n'; // [6,9]

    // kmer_hash with gaps, hashes: [2,7,8,14,1], minimiser: [2,1]
    seqan3::debug_stream << (text | seqan3::views::kmer_hash(0b101_shape) | seqan3::views::minimiser(4)) << '\n';

    /* Minimiser view with two ranges
     * The second range defines the hash values from the reverse complement, the second reverse is necessary to put the
     * hash values in the correct order. For the example here:
     * ACGTAGC | seqan3::views::complement                      => TGCATCG
     *         | std::views::reverse                            => GCTACGT
     *         | seqan3::views::kmer_hash(seqan3::ungapped{3})  => [39 (for GCA), 28 (for CTA), 49 (for TAC),
     *                                                              6 (for ACG), 27 (for CGT)]
     * "GCA" is not the reverse complement from the first k-mer in "ACGTAGC", which is "ACG", but "CGT" is.
     * Therefore, a second reverse is necessary to find the smallest value between the original sequence and its
     * reverse complement.
     */
    auto reverse_complement_hashes = text | seqan3::views::complement | std::views::reverse
                                   | seqan3::views::kmer_hash(seqan3::ungapped{3}) | std::views::reverse;
    seqan3::debug_stream << reverse_complement_hashes << '\n'; // [27,6,49,28,39]

    auto minimiser2 = seqan3::detail::minimiser_view{hashes, reverse_complement_hashes, 4};
    seqan3::debug_stream << minimiser2 << '\n'; // [6,6]
}
