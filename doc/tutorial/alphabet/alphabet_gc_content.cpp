// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [exercise]
#include <array>     // std::array
#include <iostream>  // std::cerr, std::endl
#include <cstring>   // std::strlen
#include <vector>    // std::vector

#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main (int argc, char * argv[])
{
    // Test for the presence of an input sequence.
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <sequence>" << std::endl;
        return 1;
    }

    // Get the length of the input.
    size_t len = std::strlen(argv[1]);

    // Convert the input to a dna4 sequence.
    std::vector<dna4> sequence(len);
    for (size_t idx = 0ul; idx < len; ++idx)
        assign_char_strict(sequence[idx], argv[1][idx]); // unknown characters raise an error

    // Initialise an array with count values for dna4 symbols.
    std::array<size_t, dna4::value_size> count;
    count.fill(0ul);

    // Increase the symbol count according to the sequence.
    for (dna4 symbol : sequence)
        ++count[symbol.to_rank()];

    // Calculate the GC content: (#G + #C) / (#A + #T + #G + #C).
    size_t gc = count['C'_dna4.to_rank()] + count['G'_dna4.to_rank()];
    float gc_content = 1.0f * gc / len;

    debug_stream << "The GC content of " << sequence << " is " << 100 * gc_content << "%." << std::endl;

    return 0;
}
//! [exercise]
