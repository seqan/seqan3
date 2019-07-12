// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides test utilities for seqan3::simd::simd_type types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/simd/all.hpp>

//!\cond DEV
/*!\brief #SIMD_EQ checks if the sizes and the content of two given
 * seqan3::simd::simd_type variables matches. It is like  #EXPECT_EQ, but for seqan3::simd::simd_type
 * types.
 * \ingroup simd
 * \param  left  of type seqan3::simd::simd_type
 * \param  right of type seqan3::simd::simd_type
 *
 * \attention
 * This macro can handle multiple "," which is normally a limitation of macros.
 *
 * ###Example
 *
 * \include test/snippet/core/simd/simd_test_utility.cpp
 */
#define SIMD_EQ(...) do { \
    auto [left, right] = std::make_tuple(__VA_ARGS__); \
    static_assert(seqan3::simd::Simd<decltype(left)>, "The left argument of SIMD_EQ is not a simd_type"); \
    static_assert(seqan3::simd::Simd<decltype(right)>, "The right argument of SIMD_EQ is not a simd_type"); \
    static_assert(std::is_same_v<decltype(left), decltype(right)>, "The left and right argument of SIMD_EQ don't have the same type."); \
    using _simd_traits_t = seqan3::simd::simd_traits<decltype(left)>; \
    std::vector<typename _simd_traits_t::scalar_type> left_simd(_simd_traits_t::length), right_simd(_simd_traits_t::length); \
    for (size_t i = 0; i < _simd_traits_t::length; ++i) \
    { left_simd[i] = left[i]; right_simd[i] = right[i]; } \
    EXPECT_EQ(left_simd, right_simd); \
} while (false)
//!\endcond
