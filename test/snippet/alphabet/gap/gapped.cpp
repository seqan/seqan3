// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::gapped<seqan3::dna4> gapped_letter{};
    seqan3::gapped<seqan3::dna4> converted_letter{'C'_dna4};
    seqan3::gapped<seqan3::dna4> gap_letter{seqan3::gap{}};

    seqan3::gapped<seqan3::dna4>{}.assign_char('C');
    seqan3::gapped<seqan3::dna4>{}.assign_char('-'); // gap character
    seqan3::gapped<seqan3::dna4>{}.assign_char('K'); // unknown characters map to the default/unknown
                                                     // character of the given alphabet type (i.e. A of dna4)
}
