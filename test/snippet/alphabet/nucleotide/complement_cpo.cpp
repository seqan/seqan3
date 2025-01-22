// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

using namespace seqan3::literals;

int main()
{
    auto r1 = 'A'_rna5.complement();        // calls member function rna5::complement(); r1 == 'U'_rna5
    auto r2 = seqan3::complement('A'_rna5); // calls global complement() function on the rna5 object; r2 == 'U'_rna5
}
