// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{0b1011'1000'1111};

    // begin() refers to the rightmost position.
    for (auto it = t1.begin(); it != t1.end(); ++it)
        seqan3::debug_stream << *it; // 1111'0001'1101

    seqan3::debug_stream << '\n';
}
