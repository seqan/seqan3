// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

int main()
{
    std::vector<uint16_t> memory{0, 1, 2, 3, 4, 5, 6, 7};
    uint16x8_t a = seqan3::simd::load<uint16x8_t>(memory.data());
    seqan3::debug_stream << a << '\n'; // [0,1,2,3,4,5,6,7]
    return 0;
}
