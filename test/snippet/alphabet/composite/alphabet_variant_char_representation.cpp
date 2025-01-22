// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5> var;
    var.assign_char('A'); // will be in the "dna4-state"
    var = 'A'_dna5;       // will be in the "dna5-state"
}
