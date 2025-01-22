// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<int> sequences{1, 2, 3, 4, 5, 6, 7, 8};

    seqan3::debug_stream << (sequences | seqan3::views::chunk(3)) << '\n'; // [[1,2,3],[4,5,6],[7,8]]

    // Note the behaviour on pure input ranges:
    // When incrementing the chunk view iterator, the former chunk is fully consumed by the iterator.
    // In other words, it is ensured that you always start at the beginning of a new chunk.
    auto input_view = sequences | seqan3::views::single_pass_input | seqan3::views::chunk(3);

    for (auto && val : input_view)
        seqan3::debug_stream << "first value in subrange: " << *val.begin() << '\n'; // only access first value of chunk
    // prints:
    // first value in subrange: 1
    // first value in subrange: 4
    // first value in subrange: 7
}
