// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>

using int16x8_t = seqan3::simd::simd_type_t<int16_t, 8>;
using int32x4_t = seqan3::simd::simd_type_t<int32_t, 4>;

int main()
{
    int16x8_t a{0, -1, 2, -3, 4, -5, 6, -7};
    int32x4_t b = seqan3::simd::upcast<int32x4_t>(a);

    seqan3::debug_stream << b << '\n'; // [0,-1,2,-3]
    return 0;
}
