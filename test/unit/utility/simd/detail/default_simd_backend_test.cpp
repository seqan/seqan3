// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <type_traits>

#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/detail/builtin_simd.hpp>
#include <seqan3/utility/simd/detail/default_simd_backend.hpp>

TEST(default_simd_backend, test)
{
    EXPECT_TRUE(
        (std::is_same_v<seqan3::detail::default_simd_backend<int16_t, 8>, seqan3::detail::builtin_simd<int16_t, 8>>));
}
