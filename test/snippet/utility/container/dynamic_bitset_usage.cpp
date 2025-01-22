// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{"10001100"};        // Construct from string.
    seqan3::dynamic_bitset const t2{0b1011'1000}; // Construct from binary literal or integer.

    seqan3::debug_stream << t1 << '\n'; // Debug printing inserts separators: 1000'1100
    std::cout << t2 << '\n';            // Standard printing does not: 10111000

    t1 &= t2;                           // Assign t1 the result of a binary AND of t1 and t2.
    seqan3::debug_stream << t1 << '\n'; // 1000'1000
}
