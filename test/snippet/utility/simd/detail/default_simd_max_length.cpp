// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>

int main()
{
    constexpr auto max_length = seqan3::detail::default_simd_max_length<seqan3::detail::default_simd_backend>;
    using uint8_simd_t = seqan3::simd::simd_type_t<uint8_t, max_length == 0 ? 1 : max_length>;

    uint8_simd_t a = seqan3::fill<uint8_simd_t>(4);
    uint8_simd_t b = seqan3::fill<uint8_simd_t>(5);
    uint8_simd_t c = a + b;

    seqan3::debug_stream << c << "\n"; // [9,9,9,9,9,9,9,9]
    return 0;
}
