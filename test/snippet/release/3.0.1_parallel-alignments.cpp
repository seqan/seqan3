// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alignment/pairwise/all.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using seqan3::operator""_dna4;

    auto sequence1{"ACCA"_dna4};
    auto sequence2{"ATTA"_dna4};

    seqan3::configuration alignment_config =
        seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::parallel{4};

    for (auto const & res : seqan3::align_pairwise(std::tie(sequence1, sequence2), alignment_config))
        std::cout << "Score: " << res.score() << '\n';
}
