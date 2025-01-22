// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/zip.hpp>

int main()
{
    using namespace seqan3::literals;

    // You can score two letters:
    seqan3::nucleotide_scoring_scheme scheme; // hamming is default
    seqan3::debug_stream << "Score between DNA5 A and G: " << (int)scheme.score('A'_dna5, 'G'_dna5) << "\n"; // == -1
    seqan3::debug_stream << "Score between DNA5 A and A: " << (int)scheme.score('A'_dna5, 'A'_dna5) << "\n"; // == 0

    // You can also score differenct nucleotides:
    scheme.set_simple_scheme(seqan3::match_score{3}, seqan3::mismatch_score{-2});
    seqan3::debug_stream << "Score between DNA5 A and RNA15 G: " << (int)scheme.score('A'_dna5, 'G'_rna15)
                         << "\n"; // == -2
    seqan3::debug_stream << "Score between DNA5 A and RNA15 A: " << (int)scheme.score('A'_dna5, 'A'_rna15)
                         << "\n"; // == 3

    // You can "edit" a given matrix directly:
    seqan3::nucleotide_scoring_scheme scheme2; // hamming distance is default
    seqan3::debug_stream << "Score between DNA A and G before edit: " << (int)scheme2.score('A'_dna15, 'G'_dna15)
                         << "\n"; // == -1
    scheme2.score('A'_dna15, 'G'_dna15) = 3;
    seqan3::debug_stream << "Score after editing: " << (int)scheme2.score('A'_dna15, 'G'_dna15) << "\n"; // == 3

    // You can score two sequences:
    std::vector<seqan3::dna15> one = "AGAATA"_dna15;
    std::vector<seqan3::dna15> two = "ATACTA"_dna15;
    seqan3::nucleotide_scoring_scheme scheme3; // hamming distance is default

    int score = 0;
    for (auto pair : seqan3::views::zip(one, two))
        score += scheme3.score(std::get<0>(pair), std::get<1>(pair));
    seqan3::debug_stream << "Score: " << score << "\n"; // == 0 - 1 + 0 - 1 + 0 + 0 = -2
}
