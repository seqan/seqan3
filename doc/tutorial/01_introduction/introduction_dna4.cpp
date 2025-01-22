// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//! [debug]
#include <iostream> // for std::cerr
#include <vector>   // for std::vector

#include <seqan3/alphabet/all.hpp>             // for all alphabet-related stuff
#include <seqan3/alphabet/nucleotide/dna4.hpp> // for only dna4
#include <seqan3/core/debug_stream.hpp>        // for debug_stream

int main()
{
    std::vector<seqan3::dna4> vec{seqan3::assign_char_to('A', seqan3::dna4{}),
                                  seqan3::assign_char_to('C', seqan3::dna4{}),
                                  seqan3::assign_char_to('G', seqan3::dna4{}),
                                  seqan3::assign_char_to('T', seqan3::dna4{})};
    // The above is a little cumbersome because we don't allow implicit conversions between our alphabets and `char`.
    // There is a more convenient way:
    using namespace seqan3::literals; // Lets you use operator ""_dna4 among others
    auto vec2 = "ACGT"_dna4;

    seqan3::debug_stream << vec << '\n';  // => ACGT
    seqan3::debug_stream << vec2 << '\n'; // => ACGT
}
//! [debug]
