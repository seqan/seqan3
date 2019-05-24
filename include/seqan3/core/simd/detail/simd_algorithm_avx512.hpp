// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithm implementation for SSE4.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <immintrin.h>

#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>

namespace seqan3::detail
{

/*!\brief Implementation of seqan3::simd::unpack_hi for avx512.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[in] first   The vector whose values come before the `second`.
 * \param[in] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 64
//!\endcond
inline simd_t unpack_hi(simd_t const & first, simd_t const & second)
{
    constexpr size_t scalar_size = sizeof(typename simd_traits<simd_t>::scalar_type);
    if constexpr (scalar_size == 1)
    {
        //   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
        // 63b, 63a, 62b, 62a, 62b, 61a, 60b, 60a, 59b, 59a, 58b, 58a, 57b, 57a, 56b, 56a,
        // 55b, 55a, 54b, 54a, 53b, 53a, 52b, 52a, 51b, 51a, 50b, 50a, 49b, 49a, 48b, 48a,
        // 47b, 47a, 46b, 46a, 45b, 45a, 44b, 44a, 43b, 43a, 42b, 42a, 41b, 41a, 40b, 40a,
        // 39b, 39a, 38b, 38a, 37b, 37a, 36b, 36a, 35b, 35a, 34b, 34a, 33b, 33a, 32b, 32a
        #if !defined(__clang__) && defined(__GNUC__) && __GNUC__ < 9
        __m512i idx = reinterpret_cast<__m512i>(simd_t
        {
          //  0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
            0x20, 0x60, 0x21, 0x61, 0x22, 0x62, 0x23, 0x63, 0x23, 0x63, 0x24, 0x64, 0x25, 0x65, 0x26, 0x66,
            0x28, 0x68, 0x29, 0x69, 0x2a, 0x6a, 0x2b, 0x6b, 0x2c, 0x6c, 0x2d, 0x6d, 0x2e, 0x6e, 0x2f, 0x6f,
            0x30, 0x70, 0x31, 0x71, 0x32, 0x72, 0x33, 0x73, 0x34, 0x74, 0x35, 0x75, 0x36, 0x76, 0x37, 0x77,
            0x38, 0x78, 0x39, 0x79, 0x3a, 0x7a, 0x3b, 0x7b, 0x3c, 0x7c, 0x3d, 0x75, 0x3e, 0x7e, 0x3f, 0x7f
        });
        #else // ^^^ workaround / no workaround vvv
        __m512i idx = _mm512_set_epi8
        (
            0x7f, 0x3f, 0x7e, 0x3e, 0x75, 0x3d, 0x7c, 0x3c, 0x7b, 0x3b, 0x7a, 0x3a, 0x79, 0x39, 0x78, 0x38,
            0x77, 0x37, 0x76, 0x36, 0x75, 0x35, 0x74, 0x34, 0x73, 0x33, 0x72, 0x32, 0x71, 0x31, 0x70, 0x30,
            0x6f, 0x2f, 0x6e, 0x2e, 0x6d, 0x2d, 0x6c, 0x2c, 0x6b, 0x2b, 0x6a, 0x2a, 0x69, 0x29, 0x68, 0x28,
            0x66, 0x26, 0x65, 0x25, 0x64, 0x24, 0x63, 0x23, 0x63, 0x23, 0x62, 0x22, 0x61, 0x21, 0x60, 0x20
        );
        #endif // !defined(__clang__) && defined(__GNUC__) && __GNUC__ < 9

        // TODO: does not work on clang and gcc!
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi8(reinterpret_cast<__m512i const &>(first),
                                                                 idx,
                                                                 reinterpret_cast<__m512i const &>(second)));
    }
    else if constexpr (scalar_size == 2)
    {
        //   0    1    2    3    4    5    6    7
        // 31b, 31a, 30b, 30a, 29b, 29a, 28b, 28a,
        // 27b, 27a, 26b, 26a, 25b, 25a, 24b, 24a,
        // 23b, 23a, 22b, 22a, 21b, 21a, 20b, 20a,
        // 19b, 19a, 18b, 18a, 17b, 17a, 16b, 16a
        #if !defined(__clang__) && defined(__GNUC__) && __GNUC__ < 9
        __m512i idx = reinterpret_cast<__m512i>(simd_t
        {  //  0     1     2     3     4     5     6     7
            0x10, 0x30, 0x11, 0x31, 0x12, 0x32, 0x13, 0x33,
            0x14, 0x34, 0x15, 0x35, 0x16, 0x36, 0x17, 0x37,
            0x18, 0x38, 0x19, 0x39, 0x1a, 0x3a, 0x1b, 0x3b,
            0x1c, 0x3c, 0x1d, 0x3d, 0x1e, 0x3e, 0x1f, 0x3f
        });
        #else // ^^^ workaround / no workaround vvv
        __m512i idx = _mm512_set_epi16
        (
            0x3f, 0x1f, 0x3e, 0x1e, 0x3d, 0x1d, 0x3c, 0x1c,
            0x3b, 0x1b, 0x3a, 0x1a, 0x39, 0x19, 0x38, 0x18,
            0x37, 0x17, 0x36, 0x16, 0x35, 0x15, 0x34, 0x14,
            0x33, 0x13, 0x32, 0x12, 0x31, 0x11, 0x30, 0x10
        );
        #endif // !defined(__clang__) && defined(__GNUC__) && __GNUC__ < 9
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi16(reinterpret_cast<__m512i const &>(first),
                                                                  idx,
                                                                  reinterpret_cast<__m512i const &>(second)));
    }
    else if constexpr (scalar_size == 4)
    {
        //   0    1    2    3    4    5    6    7
        // 15b, 15a, 14b, 14a, 13b, 13a, 12b, 12a,
        // 11b, 11a, 10b, 10a,  9b,  9a,  8b,  8a
        __m512i idx = _mm512_set_epi32
        (  //  0     1     2     3     4     5     6     7
            0x1f, 0x0f, 0x1e, 0x0e, 0x1d, 0x0d, 0x1c, 0x0c,
            0x1b, 0x0b, 0x1a, 0x0a, 0x19, 0x09, 0x18, 0x08
        );
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi32(reinterpret_cast<__m512i const &>(first),
                                                                  idx,
                                                                  reinterpret_cast<__m512i const &>(second)));
    }
    else if constexpr (scalar_size == 8)
    {
        //  0   1   2   3
        // 7b, 7a, 6b, 6a,
        // 5b, 5a, 4b, 4a
        __m512i idx = _mm512_set_epi64
        (  //  0     1     2     3
            0x0f, 0x07, 0x0e, 0x06,
            0x0d, 0x05, 0x0c, 0x04
        );
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi64(reinterpret_cast<__m512i const &>(first),
                                                                  idx,
                                                                  reinterpret_cast<__m512i const &>(second)));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }
}

/*!\brief Implementation of seqan3::simd::unpack_lo for avx512.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[in] first   The vector whose values come before the `second`.
 * \param[in] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 64
//!\endcond
constexpr simd_t unpack_lo(simd_t const & first, simd_t const & second)
{
    constexpr size_t scalar_size = sizeof(typename simd_traits<simd_t>::scalar_type);

    if constexpr (scalar_size == 1)
    {
        //   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
        // 31b, 31a, 30b, 30a, 29b, 29a, 28b, 28a, 27b, 27a, 26b, 26a, 25b, 25a, 24b, 24a,
        // 23b, 23a, 22b, 22a, 21b, 21a, 20b, 20a, 19b, 19a, 18b, 18a, 17b, 17a, 16b, 16a,
        // 15b, 15a, 14b, 14a, 13b, 13a, 12b, 12a, 11b, 11a, 10b, 10a, 09b, 09a, 08b, 08a,
        // 07b, 07a, 06b, 06a, 05b, 05a, 04b, 04a, 03b, 03a, 02b, 02a, 01b, 01a, 00b, 00a
        __m512i idx = reinterpret_cast<__m512i>(simd_t
        {  //  0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15
            0x5f, 0x1f, 0x5e, 0x1e, 0x5d, 0x1d, 0x5c, 0x1c, 0x5b, 0x1b, 0x5a, 0x1a, 0x59, 0x19, 0x58, 0x18,
            0x57, 0x17, 0x56, 0x16, 0x55, 0x15, 0x54, 0x14, 0x53, 0x13, 0x52, 0x12, 0x51, 0x11, 0x50, 0x10,
            0x4f, 0x0f, 0x4e, 0x0e, 0x4d, 0x0d, 0x4c, 0x0c, 0x4b, 0x0b, 0x4a, 0x0a, 0x49, 0x09, 0x48, 0x08,
            0x46, 0x06, 0x45, 0x05, 0x44, 0x04, 0x43, 0x03, 0x43, 0x03, 0x42, 0x02, 0x41, 0x01, 0x40, 0x00
        });
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi8(reinterpret_cast<__m512i const &>(lhs),
                                                                 idx,
                                                                 reinterpret_cast<__m512i const &>(second)));
    }
    else if constexpr (scalar_size == 2)
    {
        //   0    1    2    3    4    5    6    7
        // 15b, 15a, 14b, 14a, 13b, 13a, 12b, 12a,
        // 11b, 11a, 10b, 10a, 09b, 09a, 08b, 08a,
        // 07b, 07a, 06b, 06a, 05b, 05a, 04b, 04a,
        // 03b, 03a, 02b, 02a, 01b, 01a, 00b, 00a
        __m512i idx = reinterpret_cast<__m512i>(simd_t
        {  //  0     1     2     3     4     5     6     7
            0x2f, 0x0f, 0x2e, 0x0e, 0x2d, 0x0d, 0x2c, 0x0c,
            0x2b, 0x0b, 0x2a, 0x0a, 0x29, 0x09, 0x28, 0x08,
            0x27, 0x07, 0x26, 0x06, 0x25, 0x05, 0x24, 0x04,
            0x23, 0x03, 0x22, 0x02, 0x21, 0x01, 0x20, 0x00
        });
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi16(reinterpret_cast<__m512i const &>(lhs),
                                                                  idx,
                                                                  reinterpret_cast<__m512i const &>(second)));
    }
    else if constexpr (scalar_size == 4)
    {
        //   0  1   2   3   4   5   6   7
        // 7b, 7a, 6b, 6a, 5b, 5a, 4b, 4a,
        // 3b, 3a, 2b, 2a, 1b, 1a, 0b, 0a
        __m512i idx = _mm512_set_epi32
        (  //  0     1     2     3     4     5     6     7
            0x17, 0x07, 0x16, 0x06, 0x15, 0x05, 0x14, 0x04,
            0x13, 0x03, 0x12, 0x02, 0x11, 0x01, 0x10, 0x00
        );
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi32(reinterpret_cast<__m512i const &>(first),
                                                                  idx,
                                                                  reinterpret_cast<__m512i const &>(second)));
    }
    else if constexpr (scalar_size == 8)
    {
        //  0   1   2   3
        // 3b, 3a, 2b, 2a,
        // 1b, 1a, 0b, 0a
        __m512i idx = _mm512_set_epi64
        (  //  0     1     2     3
            0x0b, 0x03, 0x0a, 0x02,
            0x09, 0x01, 0x08, 0x00
        );
        return reinterpret_cast<simd_t>(_mm512_permutex2var_epi64(reinterpret_cast<__m512i const &>(first),
                                                                  idx,
                                                                  reinterpret_cast<__m512i const &>(second)));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }
}

} // namespace seqan3::detail
