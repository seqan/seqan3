// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{0b1011'1000'1111};

    // Positions are 0-based and start at the rightmost bit.
    for (size_t i = 0; i < t1.size(); ++i)
        seqan3::debug_stream << t1[i]; // 1111'0001'1101

    seqan3::debug_stream << '\n';
}
