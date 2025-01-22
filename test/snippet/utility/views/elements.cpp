// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/aliases.hpp> // includes seqan3::dna4q
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/elements.hpp>

int main()
{
    using namespace seqan3::literals;

    // Create a vector of dna4 quality composite alphabet.
    std::vector<seqan3::dna4q> qv{{'A'_dna4, '0'_phred42},
                                  {'C'_dna4, '1'_phred42},
                                  {'G'_dna4, '2'_phred42},
                                  {'T'_dna4, '3'_phred42}};

    seqan3::debug_stream << (qv | seqan3::views::elements<0> | seqan3::views::to_char) << '\n'; // Prints [A,C,G,T]
    seqan3::debug_stream << (qv | seqan3::views::elements<1> | seqan3::views::to_char) << '\n'; // Prints [!,",#,$]
}
