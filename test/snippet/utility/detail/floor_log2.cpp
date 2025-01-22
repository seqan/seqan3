// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/math.hpp>

int main()
{
    // the same as std::floor(std::log2(x)), but exact for unsigned integers
    seqan3::debug_stream << "floor_log2(2^0 + 0): " << seqan3::detail::floor_log2(1u) << '\n'; // 0u
    seqan3::debug_stream << "floor_log2(2^1 + 0): " << seqan3::detail::floor_log2(2u) << '\n'; // 1u
    seqan3::debug_stream << "floor_log2(2^1 + 1): " << seqan3::detail::floor_log2(3u) << '\n'; // 1u
    seqan3::debug_stream << "floor_log2(2^2 + 0): " << seqan3::detail::floor_log2(4u) << '\n'; // 2u
    seqan3::debug_stream << "floor_log2(2^2 + 1): " << seqan3::detail::floor_log2(5u) << '\n'; // 2u
    seqan3::debug_stream << "floor_log2(2^2 + 2): " << seqan3::detail::floor_log2(6u) << '\n'; // 2u
    seqan3::debug_stream << "floor_log2(2^2 + 3): " << seqan3::detail::floor_log2(7u) << '\n'; // 2u
    seqan3::debug_stream << "floor_log2(2^3 + 0): " << seqan3::detail::floor_log2(8u) << '\n'; // 3u
    seqan3::debug_stream << "floor_log2(2^3 + 1): " << seqan3::detail::floor_log2(9u) << '\n'; // 3u

    return 0;
}
