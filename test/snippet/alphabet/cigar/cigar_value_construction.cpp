// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    // creates 10M, as the cigar_op field is not provided.
    seqan3::cigar letter1{10};
    seqan3::debug_stream << "letter1: " << letter1 << '\n'; // 10M

    // creates 0I, as the integer field is not provided.
    seqan3::cigar letter2{'I'_cigar_operation};
    seqan3::debug_stream << "letter2: " << letter2 << '\n'; // 0I

    // creates 10I, as both fields are explicitly given.
    seqan3::cigar letter3{10, 'I'_cigar_operation};
    seqan3::debug_stream << "letter3: " << letter3 << '\n'; // 10I
}
