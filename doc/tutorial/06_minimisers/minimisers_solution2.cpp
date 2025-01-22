// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};

    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 4 is used. The seed is set
    // to 0, so lexicographical ordering is used.
    auto example_a =
        text
        | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{4}, seqan3::seed{0});
    // results in: [81, 70, 27, 109, 97, 216, 97, 109, 26, 22, 5]
    // representing the k-mers [CCAC, CACG, ACGT, CGTC, cgac, TCGA, CGAC, cgtc, ACGG, accg, aacc]
    seqan3::debug_stream << example_a << '\n';

    auto example_b =
        text
        | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{8}, seqan3::seed{0});
    // results in: [27, 97, 26, 22, 5] representing the k-mers [ACGT, CGAC, ACGG, accg, aacc]
    seqan3::debug_stream << example_b << '\n';

    auto example_c = text | seqan3::views::minimiser_hash(0b1'0101_shape, seqan3::window_size{8}, seqan3::seed{0});
    // results in: [9, 18, 7, 6] representing the k-mers [A.G.C, C.A.G, a.c.t, a.c.g]
    seqan3::debug_stream << example_c << '\n';
}
