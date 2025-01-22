// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

int main()
{
    uint16x8_t a = seqan3::fill<uint16x8_t>(4);
    seqan3::debug_stream << a << "\n"; // [4,4,4,4,4,4,4,4]

    // or:

    uint16x8_t b = seqan3::simd::fill<uint16x8_t>(4);
    seqan3::debug_stream << b << "\n"; // [4,4,4,4,4,4,4,4]
    return 0;
}
