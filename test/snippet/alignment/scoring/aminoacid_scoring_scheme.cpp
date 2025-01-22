// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/zip.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aminoacid_scoring_scheme scheme{seqan3::aminoacid_similarity_matrix::blosum62};
    // How to score two letters:
    seqan3::debug_stream << "blosum62 score for T and S: " << (int)scheme.score('T'_aa27, 'S'_aa27) << "\n"; // == 1

    scheme.set_similarity_matrix(seqan3::aminoacid_similarity_matrix::blosum80);
    // You can also score aa20 against aa27:
    seqan3::debug_stream << "blosum80 score for 'T'_aa27 and 'S'_aa20: " << (int)scheme.score('T'_aa27, 'S'_aa20)
                         << "\n"; // == 2
    scheme.set_hamming_distance();
    seqan3::debug_stream << "Hamming distance between T and S: " << (int)scheme.score('T'_aa27, 'S'_aa20)
                         << "\n"; // == -1
    seqan3::debug_stream << "Hamming distance between T and T: " << (int)scheme.score('T'_aa27, 'T'_aa20)
                         << "\n"; // == 0

    seqan3::aminoacid_scoring_scheme scheme2{seqan3::aminoacid_similarity_matrix::blosum80};
    // You can "edit" a given matrix directly:
    seqan3::debug_stream << "blosum80 score between T and S: " << (int)scheme2.score('T'_aa27, 'S'_aa27)
                         << "\n"; // == 2
    auto & cell = scheme2.score('T'_aa27, 'S'_aa27);
    cell = 3;
    seqan3::debug_stream << "New score after editing entry: " << (int)scheme2.score('T'_aa27, 'S'_aa27) << "\n"; // == 3

    std::vector<seqan3::aa27> one = "ALIGATOR"_aa27;
    std::vector<seqan3::aa27> two = "ANIMATOR"_aa27;

    seqan3::aminoacid_scoring_scheme scheme3{seqan3::aminoacid_similarity_matrix::blosum62};
    // You can also score two sequences:
    int score = 0;
    for (auto pair : seqan3::views::zip(one, two))
        score += scheme3.score(std::get<0>(pair), std::get<1>(pair));
    seqan3::debug_stream << "Score: " << score << "\n"; // 4 + -3 + 4 + -3 + 4 + 5 + -1 + 5 = 15
}
