// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides test utilities for seqan3::simd::simd_type types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <array>

#include <seqan3/utility/simd/concept.hpp>

//!\cond DEV
/*!\brief #SIMD_EQ checks if the sizes and the content of two given
 * seqan3::simd::simd_type variables matches. It is like  #EXPECT_EQ, but for seqan3::simd::simd_type
 * types.
 * \ingroup utility_simd
 * \param[in] left_simd_argument Left hand side operand of type seqan3::simd::simd_type.
 * \param[in] right_simd_argument Right hand side operand of type seqan3::simd::simd_type.
 *
 * \attention
 * This macro can handle multiple "," which is normally a limitation of macros.
 *
 * ###Example
 *
 * \include test/snippet/test/simd_utility.cpp
 */
#define SIMD_EQ(...)                                                                                                   \
    do                                                                                                                 \
    {                                                                                                                  \
        auto [left_simd_argument, right_simd_argument] = std::tuple{__VA_ARGS__};                                      \
        static_assert(seqan3::simd::simd_concept<decltype(left_simd_argument)>,                                        \
                      "The left argument of SIMD_EQ is not a simd_type");                                              \
        static_assert(seqan3::simd::simd_concept<decltype(right_simd_argument)>,                                       \
                      "The right argument of SIMD_EQ is not a simd_type");                                             \
        static_assert(std::is_same_v<decltype(left_simd_argument), decltype(right_simd_argument)>,                     \
                      "The left and right argument of SIMD_EQ don't have the same type.");                             \
        using _simd_traits_t = seqan3::simd::simd_traits<decltype(left_simd_argument)>;                                \
        std::array<typename _simd_traits_t::scalar_type, _simd_traits_t::length> left_simd_argument_as_array{},        \
            right_simd_argument_as_array{};                                                                            \
        for (size_t i = 0; i < _simd_traits_t::length; ++i)                                                            \
        {                                                                                                              \
            left_simd_argument_as_array[i] = left_simd_argument[i];                                                    \
            right_simd_argument_as_array[i] = right_simd_argument[i];                                                  \
        }                                                                                                              \
        EXPECT_EQ(left_simd_argument_as_array, right_simd_argument_as_array);                                          \
    }                                                                                                                  \
    while (false)
//!\endcond
