// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>

// The simd type with 8 unsigned shorts.
using uint16x8_t = seqan3::simd_type_t<uint16_t, 8>;

int main()
{
    // Generate ascending simd index using a iota view in combination with a transform view.
    auto simd_iota_view = std::views::iota(0, 10)
                        | std::views::transform(
                              [](uint16_t const idx)
                              {
                                  return seqan3::simd::fill<uint16x8_t>(idx);
                              });

    for (auto && simd_id : simd_iota_view)
        seqan3::debug_stream << simd_id << '\n'; // [0, 0, ..., 0], [1, 1, ..., 1], ... [9, 9, ..., 9]

    return 0;
}
