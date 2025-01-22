// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges> // For std::ranges::size()

#include <seqan3/core/debug_stream.hpp>       // For seqan3::debug_stream
#include <seqan3/search/kmer_index/shape.hpp> // For seqan3::shape

using namespace seqan3::literals;

int main()
{
    seqan3::shape s0{seqan3::ungapped{5}};                 // represents "11111", i.e. ungapped 5-mer
    seqan3::debug_stream << s0 << '\n';                    // prints "[1,1,1,1,1]"
    seqan3::debug_stream << std::ranges::size(s0) << '\n'; // prints "5"

    seqan3::shape s1{seqan3::bin_literal{0b101}};          // represents "101", i.e. gapped 3-mer
    seqan3::debug_stream << s1 << '\n';                    // prints "[1,0,1]"
    seqan3::debug_stream << std::ranges::size(s1) << '\n'; // prints "3"

    seqan3::shape s2{0b101_shape};                         // Same as previous
    seqan3::debug_stream << s2 << '\n';                    // prints "[1,0,1]"
    seqan3::debug_stream << std::ranges::size(s2) << '\n'; // prints "3"
}
