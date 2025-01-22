// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>

// The simd type with 8 unsigned shorts.
using uint16x8_t = seqan3::simd_type_t<uint16_t, 8>;

int main()
{
    // Generate ascending simd index.
    for (auto && simd_id : seqan3::views::iota_simd<uint16x8_t>(0, 10))
        seqan3::debug_stream << simd_id << '\n'; // [0, 0, ..., 0], [1, 1, ..., 1], ... [9, 9, ..., 9]

    return 0;
}
