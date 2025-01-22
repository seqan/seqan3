// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <utility>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/hamming_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector s1 = "ACGTGAACTGACT"_dna4;
    seqan3::dna4_vector s2 = "ACGAAGACCGAT"_dna4;

    // Configure the alignment kernel.
    auto config =
        seqan3::align_cfg::method_global{} | seqan3::align_cfg::scoring_scheme{seqan3::hamming_scoring_scheme{}};

    // Invoke the pairwise alignment which returns a lazy range over alignment results.
    auto results = seqan3::align_pairwise(std::tie(s1, s2), config);
    auto & res = *results.begin();
    seqan3::debug_stream << "Score: " << res.score() << '\n';
}
