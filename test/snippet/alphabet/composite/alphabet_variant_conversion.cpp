// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::alphabet_variant<seqan3::dna4, seqan3::gap> letter1{'C'_rna4};
    seqan3::alphabet_variant<seqan3::dna4, seqan3::gap> letter2 = 'C'_rna4;
}
