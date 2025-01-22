// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/simd/all.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

static_assert(std::is_same_v<seqan3::simd::simd_traits<uint16x8_t>::scalar_type, uint16_t>);
static_assert(seqan3::simd::simd_traits<uint16x8_t>::length == 8);
static_assert(seqan3::simd::simd_traits<uint16x8_t>::max_length == 16);
