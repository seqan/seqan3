// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges> // For std::ranges::size()

#include <seqan3/core/debug_stream.hpp>       // For seqan3::debug_stream
#include <seqan3/search/kmer_index/shape.hpp> // For seqan3::shape

using namespace seqan3::literals;

#if 1
// this case is fun, as seqan3::shape has no explicit std::cout/debug_stream overload.
// * if we have
//   printer_order<debug_stream_printer, std_printer, input_range_printer>
//   the std::cout overload of seqan3::dynamic_bitset will win as it is the most specific
//   seqan3::debug_stream << s0; will print 11111
// * if we have
//   printer_order<debug_stream_printer, input_range_printer, std_printer>
//   the input_range_printer will win since seqan3::dynamic_bitset is an input_range
//   seqan3::debug_stream << s0; will print [1,1,1,1,1]
//
// Interestingly seqan3::dynamic_bitset has also an debug_stream overload, but this one will
// only be active if the type is seqan3::dynamic_bitset
// seqan3::debug_stream << (seqan3::dynamic_bitset<58>&)s0; will print 1'1111
#endif

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
