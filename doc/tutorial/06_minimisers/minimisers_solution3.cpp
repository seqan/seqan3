// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};

    // This would lead to an static assert error, because the shape.size() equals the window size. (Remember, the input
    // parameter for the minimiser view is calculated by: window size - k-mer size + 1, here: 4 - 4 + 1 = 1.)
    // Therefore, kmer_hash needs to be used.
    /* auto example_a = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{4}})
                             | seqan3::views::minimiser(1); */

    // results in: [81, 70, 27, 109, 182, 216, 97, 134, 26, 107, 175]
    // representing the k-mers [CCAC, CACG, ACGT, CGTC, GTCG, TCGA, CGAC, GACG, ACGG, CGGT, GGTT]
    auto example_a = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{4}});
    seqan3::debug_stream << example_a << '\n';

    auto example_b = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{4}}) | seqan3::views::minimiser(5);
    // results in: [27, 97, 26] representing the k-mers [ACGT, CGAC, ACGG]
    seqan3::debug_stream << example_b << '\n';

    auto example_c = text | seqan3::views::kmer_hash(0b1'0101_shape) | seqan3::views::minimiser(4);
    // results in: [9, 18, 11] representing the k-mers [A.G.C, C.A.G, A.G.T]
    seqan3::debug_stream << example_c << '\n';
}
