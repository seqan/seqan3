// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include <seqan3/test/simd_utility.hpp>
#include <seqan3/utility/simd/all.hpp>

TEST(simd, simd_eq)
{
    using int16x8_t = seqan3::simd::simd_type_t<int16_t, 8u>;
    int16x8_t a{4, 3, 2, 1, 0, -1, -2, -3};

    SIMD_EQ(a, int16x8_t{4, 3, 2, 1, 0, -1, -2, -3});
}
