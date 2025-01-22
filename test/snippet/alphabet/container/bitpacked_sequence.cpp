// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> v0{"ACGT"_dna4};                // data occupies 4 bytes in memory
    seqan3::bitpacked_sequence<seqan3::dna4> v1{"ACGT"_dna4}; // data occupies 1 byte in memory
}
