// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4 letter{'C'_dna4};

    letter.assign_char('F');                // Characters other than IUPAC characters are implicitly converted to A.
    seqan3::debug_stream << letter << '\n'; // prints "A"

    // IUPAC characters are implicitly converted to their best fitting representative
    seqan3::debug_stream << letter.assign_char('R') << '\n'; // prints "A"
    seqan3::debug_stream << letter.assign_char('Y') << '\n'; // prints "C"
    seqan3::debug_stream << letter.assign_char('S') << '\n'; // prints "C"
    seqan3::debug_stream << letter.assign_char('W') << '\n'; // prints "A"
    seqan3::debug_stream << letter.assign_char('K') << '\n'; // prints "G"
    seqan3::debug_stream << letter.assign_char('M') << '\n'; // prints "A"
    seqan3::debug_stream << letter.assign_char('B') << '\n'; // prints "C"
    seqan3::debug_stream << letter.assign_char('D') << '\n'; // prints "A"
    seqan3::debug_stream << letter.assign_char('H') << '\n'; // prints "A"
    seqan3::debug_stream << letter.assign_char('V') << '\n'; // prints "A"

    letter.assign_char('a');                // Lower case letters are the same as their upper case equivalent.
    seqan3::debug_stream << letter << '\n'; // prints "A"
}
