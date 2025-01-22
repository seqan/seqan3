// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{"101110001111"}; // Use character literal.
    seqan3::dynamic_bitset t2{"000010001111"}; // Leading zeros are kept.

    seqan3::debug_stream << t1 << '\n'; // 1011'1000'1111
    seqan3::debug_stream << t2 << '\n'; // 0000'1000'1111
}
