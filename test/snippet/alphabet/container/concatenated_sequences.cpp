// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::concatenated_sequences<seqan3::dna4_vector> concat1{"ACGT"_dna4, "GAGGA"_dna4};
    seqan3::debug_stream << concat1[0] << '\n'; // "ACGT"

    std::vector<seqan3::dna4_vector> concat2{"ACTA"_dna4, "AGGA"_dna4};

    concat1 = concat2; // you can assign from other ranges

    concat2[0] = "ATTA"_dna4;                   // this works for vector of vector
    concat1[0][1] = 'T'_dna4;                   // and this works for concatenated_sequences
    seqan3::debug_stream << concat1[0] << '\n'; // "ATTA"

    // if you know that you will be adding ten vectors of length ten:
    std::vector<seqan3::dna4> vector_of_length10{"ACGTACGTAC"_dna4};

    concat1.reserve(10);
    concat1.concat_reserve(10 * vector_of_length10.size());
    while (concat1.size() < 10)
    {
        // ...
        concat1.push_back(vector_of_length10);
    }
}
