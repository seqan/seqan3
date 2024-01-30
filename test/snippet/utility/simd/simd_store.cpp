// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

int main()
{
    uint16x8_t a = seqan3::simd::iota<uint16x8_t>(0);
    std::vector<uint16_t> memory(seqan3::simd::simd_traits<uint16x8_t>::length);
    seqan3::simd::store(memory.data(), a);
    seqan3::debug_stream << memory << '\n'; // [0,1,2,3,4,5,6,7]
    return 0;
}
