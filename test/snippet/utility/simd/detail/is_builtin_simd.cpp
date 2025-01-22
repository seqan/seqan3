// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/simd/detail/builtin_simd.hpp>

using int8x16_t = seqan3::detail::builtin_simd<int8_t, 16>::type;

static_assert(seqan3::detail::is_builtin_simd<int8x16_t>::value);
static_assert(!seqan3::detail::is_builtin_simd<int8_t>::value);
