// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithm implementation for AVX2.
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

/*!\brief Implementation of seqan3::simd::unpack_hi for avx2.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[in] first   The vector whose values come before the `second`.
 * \param[in] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 32
//!\endcond
inline simd_t unpack_hi(simd_t const & first, simd_t const & second)
{
    constexpr size_t scalar_size = sizeof(typename simd_traits<simd_t>::scalar_type);

    [[maybe_unused]] __m256i tmp_lo{};
    [[maybe_unused]] __m256i tmp_hi{};

    if constexpr (scalar_size == 1)
    {
        tmp_lo = _mm256_unpacklo_epi8(reinterpret_cast<__m256i const &>(first),
                                      reinterpret_cast<__m256i const &>(second));
        tmp_hi = _mm256_unpackhi_epi8(reinterpret_cast<__m256i const &>(first),
                                      reinterpret_cast<__m256i const &>(second));
    }
    else if constexpr (scalar_size == 2)
    {
        tmp_lo = _mm256_unpacklo_epi16(reinterpret_cast<__m256i const &>(first),
                                       reinterpret_cast<__m256i const &>(second));
        tmp_hi = _mm256_unpackhi_epi16(reinterpret_cast<__m256i const &>(first),
                                       reinterpret_cast<__m256i const &>(second));
    }
    else if constexpr (scalar_size == 4)
    {
        tmp_lo = _mm256_unpacklo_epi32(reinterpret_cast<__m256i const &>(first),
                                       reinterpret_cast<__m256i const &>(second));
        tmp_hi = _mm256_unpackhi_epi32(reinterpret_cast<__m256i const &>(first),
                                       reinterpret_cast<__m256i const &>(second));
    }
    else if constexpr (scalar_size == 8)
    {
        tmp_lo = _mm256_unpacklo_epi64(reinterpret_cast<__m256i const &>(first),
                                       reinterpret_cast<__m256i const &>(second));
        tmp_hi = _mm256_unpackhi_epi64(reinterpret_cast<__m256i const &>(first),
                                       reinterpret_cast<__m256i const &>(second));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }

    return reinterpret_cast<simd_t>(_mm256_permute2f128_si256(tmp_lo, tmp_hi, 0x31));
}

} // namespace seqan3::detail
