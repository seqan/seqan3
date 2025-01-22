// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{0b1011'1000'1111}; // Use binary literal.
    seqan3::dynamic_bitset t2{0b0000'1000'1111}; // Leading zeros are stripped.
    seqan3::dynamic_bitset t3{832};              // Use a number.

    seqan3::debug_stream << t1 << '\n'; // 1011'1000'1111
    seqan3::debug_stream << t2 << '\n'; // 1000'1111
    seqan3::debug_stream << t3 << '\n'; // 1101'0000'00
}
