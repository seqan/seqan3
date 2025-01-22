// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::masked<seqan3::dna4> dna4_masked{};
    seqan3::masked<seqan3::dna4> dna4_another_masked{'A'_dna4, seqan3::mask::unmasked};
    // create a dna4 masked alphabet with an unmasked A

    dna4_masked.assign_char('a'); // assigns a masked 'A'_dna4

    if (dna4_masked.to_char() != dna4_another_masked.to_char())
    {
        seqan3::debug_stream << dna4_masked.to_char() << " is not the same as " << dna4_another_masked.to_char()
                             << "\n";
    }
}
