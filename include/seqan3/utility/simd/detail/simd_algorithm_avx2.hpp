// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides specific algorithm implementations for AVX2 instruction set.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>

#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/detail/builtin_simd.hpp>
#include <seqan3/utility/simd/detail/builtin_simd_intrinsics.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

//-----------------------------------------------------------------------------
// forward declare avx2 simd algorithms that use avx2 intrinsics
//-----------------------------------------------------------------------------

namespace seqan3::detail
{
/*!\copydoc seqan3::simd::load
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <simd::simd_concept simd_t>
constexpr simd_t load_avx2(void const * mem_addr);

/*!\copydoc seqan3::simd::store
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <simd::simd_concept simd_t>
constexpr void store_avx2(void * mem_addr, simd_t const & simd_vec);

/*!\copydoc seqan3::simd::transpose
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <simd::simd_concept simd_t>
inline void transpose_matrix_avx2(std::array<simd_t, simd_traits<simd_t>::length> & matrix);

/*!\copydoc seqan3::detail::upcast_signed
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed_avx2(source_simd_t const & src);

/*!\copydoc seqan3::detail::upcast_unsigned
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned_avx2(source_simd_t const & src);

/*!\copydoc seqan3::detail::extract_half
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_half_avx2(simd_t const & src);

/*!\copydoc seqan3::detail::extract_quarter
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_avx2(simd_t const & src);

/*!\copydoc seqan3::detail::extract_eighth
 * \attention This is the implementation for AVX2 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_avx2(simd_t const & src);

} // namespace seqan3::detail

//-----------------------------------------------------------------------------
// implementation
//-----------------------------------------------------------------------------

#ifdef __AVX2__

namespace seqan3::detail
{

template <simd::simd_concept simd_t>
constexpr simd_t load_avx2(void const * mem_addr)
{
    return reinterpret_cast<simd_t>(_mm256_loadu_si256(reinterpret_cast<__m256i const *>(mem_addr)));
}

template <simd::simd_concept simd_t>
constexpr void store_avx2(void * mem_addr, simd_t const & simd_vec)
{
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(mem_addr), reinterpret_cast<__m256i const &>(simd_vec));
}

template <simd::simd_concept simd_t>
inline void transpose_matrix_avx2(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    // emulate missing _mm256_unpacklo_epi128/_mm256_unpackhi_epi128 instructions
    auto _mm256_unpacklo_epi128 = [](__m256i const & a, __m256i const & b)
    {
        return _mm256_permute2x128_si256(a, b, 0x20);
    };

    auto _mm256_unpackhi_epi128 = [](__m256i const & a, __m256i const & b)
    {
        return _mm256_permute2x128_si256(a, b, 0x31);
    };

    // A look-up table to reverse the lowest 4 bits in order to permute the transposed rows.
    static uint8_t const bit_rev[] = {0,  8,  4,  12, 2,  10, 6,  14, 1,  9,  5,  13, 3,  11, 7,  15,
                                      16, 24, 20, 28, 18, 26, 22, 30, 17, 25, 21, 29, 19, 27, 23, 31};

    // transpose a 32x32 byte matrix
    __m256i tmp1[32];
    for (int i = 0; i < 16; ++i)
    {
        tmp1[i] = _mm256_unpacklo_epi8(reinterpret_cast<__m256i const &>(matrix[2 * i]),
                                       reinterpret_cast<__m256i const &>(matrix[2 * i + 1]));
        tmp1[i + 16] = _mm256_unpackhi_epi8(reinterpret_cast<__m256i const &>(matrix[2 * i]),
                                            reinterpret_cast<__m256i const &>(matrix[2 * i + 1]));
    }
    __m256i tmp2[32];
    for (int i = 0; i < 16; ++i)
    {
        tmp2[i] = _mm256_unpacklo_epi16(tmp1[2 * i], tmp1[2 * i + 1]);
        tmp2[i + 16] = _mm256_unpackhi_epi16(tmp1[2 * i], tmp1[2 * i + 1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        tmp1[i] = _mm256_unpacklo_epi32(tmp2[2 * i], tmp2[2 * i + 1]);
        tmp1[i + 16] = _mm256_unpackhi_epi32(tmp2[2 * i], tmp2[2 * i + 1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        tmp2[i] = _mm256_unpacklo_epi64(tmp1[2 * i], tmp1[2 * i + 1]);
        tmp2[i + 16] = _mm256_unpackhi_epi64(tmp1[2 * i], tmp1[2 * i + 1]);
    }
    for (int i = 0; i < 16; ++i)
    {
        matrix[bit_rev[i]] = reinterpret_cast<simd_t>(_mm256_unpacklo_epi128(tmp2[2 * i], tmp2[2 * i + 1]));
        matrix[bit_rev[i + 16]] = reinterpret_cast<simd_t>(_mm256_unpackhi_epi128(tmp2[2 * i], tmp2[2 * i + 1]));
    }
}

template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed_avx2(source_simd_t const & src)
{
    __m128i const & tmp = _mm256_castsi256_si128(reinterpret_cast<__m256i const &>(src));
    if constexpr (simd_traits<source_simd_t>::length == 32) // cast from epi8 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 16) // to epi16
            return reinterpret_cast<target_simd_t>(_mm256_cvtepi8_epi16(tmp));
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi32
            return reinterpret_cast<target_simd_t>(_mm256_cvtepi8_epi32(tmp));
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi64
            return reinterpret_cast<target_simd_t>(_mm256_cvtepi8_epi64(tmp));
    }
    else if constexpr (simd_traits<source_simd_t>::length == 16) // cast from epi16 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi32
            return reinterpret_cast<target_simd_t>(_mm256_cvtepi16_epi32(tmp));
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi64
            return reinterpret_cast<target_simd_t>(_mm256_cvtepi16_epi64(tmp));
    }
    else // cast from epi32 to epi64
    {
        static_assert(simd_traits<source_simd_t>::length == 8, "Expected 32 bit scalar type.");
        return reinterpret_cast<target_simd_t>(_mm256_cvtepi32_epi64(tmp));
    }
}

template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned_avx2(source_simd_t const & src)
{
    __m128i const & tmp = _mm256_castsi256_si128(reinterpret_cast<__m256i const &>(src));
    if constexpr (simd_traits<source_simd_t>::length == 32) // cast from epi8 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 16) // to epi16
            return reinterpret_cast<target_simd_t>(_mm256_cvtepu8_epi16(tmp));
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi32
            return reinterpret_cast<target_simd_t>(_mm256_cvtepu8_epi32(tmp));
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi64
            return reinterpret_cast<target_simd_t>(_mm256_cvtepu8_epi64(tmp));
    }
    else if constexpr (simd_traits<source_simd_t>::length == 16) // cast from epi16 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi32
            return reinterpret_cast<target_simd_t>(_mm256_cvtepu16_epi32(tmp));
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi64
            return reinterpret_cast<target_simd_t>(_mm256_cvtepu16_epi64(tmp));
    }
    else // cast from epi32 to epi64
    {
        static_assert(simd_traits<source_simd_t>::length == 8, "Expected 32 bit scalar type.");
        return reinterpret_cast<target_simd_t>(_mm256_cvtepu32_epi64(tmp));
    }
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_half_avx2(simd_t const & src)
{
    return reinterpret_cast<simd_t>(
        _mm256_castsi128_si256(_mm256_extracti128_si256(reinterpret_cast<__m256i const &>(src), index)));
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_avx2(simd_t const & src)
{
    return reinterpret_cast<simd_t>(
        _mm256_castsi128_si256(_mm_cvtsi64_si128(_mm256_extract_epi64(reinterpret_cast<__m256i const &>(src), index))));
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_avx2(simd_t const & src)
{
    return reinterpret_cast<simd_t>(
        _mm256_castsi128_si256(_mm_cvtsi32_si128(_mm256_extract_epi32(reinterpret_cast<__m256i const &>(src), index))));
}

} // namespace seqan3::detail

#endif // __AVX2__
