// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/pairwise/all.hpp>  // for seqan3::align_cfg and seqan3::align_pairwise
#include <seqan3/alignment/scoring/all.hpp>   // for seqan3::aminoacid_scoring_scheme and
                                              //     seqan3::aminoacid_similarity_matrix
#include <seqan3/alphabet/aminoacid/aa27.hpp> // for seqan3::operator""_aa27
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    auto seq1 = "QFSEEILSDIYCWMLQCGQERAV"_aa27;
    auto seq2 = "AFLPGWQEENKLSKIWMKDCGCLW"_aa27;

    // Configure the alignment kernel.
    auto config =
        seqan3::align_cfg::method_global{}
        | seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{
            seqan3::aminoacid_similarity_matrix::blosum62}}
        | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-9}, seqan3::align_cfg::extension_score{-2}};

    for (auto const & res : seqan3::align_pairwise(std::tie(seq1, seq2), config))
        seqan3::debug_stream << "Score: " << res.score() << '\n';
}
