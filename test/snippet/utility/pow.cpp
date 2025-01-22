// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/math.hpp>

int main()
{
    // Uses specialisation for signed integers.
    seqan3::debug_stream << seqan3::pow(2, 3u) << '\n';  // Prints 8
    seqan3::debug_stream << seqan3::pow(-2, 3u) << '\n'; // Prints -8

    // Uses specialisation for unsigned integers.
    seqan3::debug_stream << seqan3::pow(2u, 3u) << '\n'; // Prints 8

    // Uses `std::pow`.
    seqan3::debug_stream << seqan3::pow(2, 3) << '\n';   // Prints 8
    seqan3::debug_stream << seqan3::pow(2u, 3) << '\n';  // Prints 8
    seqan3::debug_stream << seqan3::pow(2.0, 3) << '\n'; // Prints 8

    // 5^25 should be 298023223876953125.
    seqan3::debug_stream << seqan3::pow(5u, 25u) << '\n';                     // Prints 298023223876953125
    seqan3::debug_stream << static_cast<uint64_t>(std::pow(5u, 25u)) << '\n'; // Prints 298023223876953152 (wrong!)
}
