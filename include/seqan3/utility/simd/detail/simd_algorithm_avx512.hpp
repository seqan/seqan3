// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides specific algorithm implementations for AVX512 instruction set.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>

#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/detail/builtin_simd.hpp>
#include <seqan3/utility/simd/detail/builtin_simd_intrinsics.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

//-----------------------------------------------------------------------------
// forward declare avx512 simd algorithms that use avx512 intrinsics
//-----------------------------------------------------------------------------

namespace seqan3::detail
{
/*!\copydoc seqan3::simd::load
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <simd::simd_concept simd_t>
constexpr simd_t load_avx512(void const * mem_addr);

/*!\copydoc seqan3::simd::store
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <simd::simd_concept simd_t>
constexpr void store_avx512(void * mem_addr, simd_t const & simd_vec);

/*!\copydoc seqan3::simd::transpose
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <simd::simd_concept simd_t>
inline void transpose_matrix_avx512(std::array<simd_t, simd_traits<simd_t>::length> & matrix);

/*!\copydoc seqan3::detail::upcast_signed
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed_avx512(source_simd_t const & src);

/*!\copydoc seqan3::detail::upcast_unsigned
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned_avx512(source_simd_t const & src);

/*!\copydoc seqan3::detail::extract_half
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_half_avx512(simd_t const & src);

/*!\copydoc seqan3::detail::extract_quarter
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_avx512(simd_t const & src);

/*!\copydoc seqan3::detail::extract_eighth
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_avx512(simd_t const & src);

} // namespace seqan3::detail

//-----------------------------------------------------------------------------
// implementation
//-----------------------------------------------------------------------------

#ifdef __AVX512F__

namespace seqan3::detail
{

template <simd::simd_concept simd_t>
constexpr simd_t load_avx512(void const * mem_addr)
{
    return reinterpret_cast<simd_t>(_mm512_loadu_si512(mem_addr));
}

template <simd::simd_concept simd_t>
constexpr void store_avx512(void * mem_addr, simd_t const & simd_vec)
{
    _mm512_storeu_si512(mem_addr, reinterpret_cast<__m512i const &>(simd_vec));
}

#    if defined(__AVX512BW__)
template <simd::simd_concept simd_t>
inline void transpose_matrix_avx512(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    // Transposing a 64x64 byte matrix in 6x64 instructions using AVX512 intrinsics.
    // Step 1: Unpack 8-bit operands.

    // Note: _mm512_unpack* operates on 128-bit lanes, i.e. the following pattern applies:
    //  * for lower half ('lo')
    //       d[127:0]   = interleave(a[63:0],    b[63:0])
    //       d[255:128] = interleave(a[191:128], b[191:128])
    //       d[383:256] = interleave(a[319:256], b[319:256])
    //       d[511:384] = interleave(a[447:384], b[447:384])
    //  * for higher half ('hi')
    //       d[127:0]   = interleave(a[127:64],  b[127:64])
    //       d[255:128] = interleave(a[255:192], b[255:192])
    //       d[383:256] = interleave(a[319:320], b[383:320])
    //       d[511:384] = interleave(a[511:448], b[511:448])
    __m512i tmp1[64];
    for (int i = 0; i < 32; ++i)
    {
        tmp1[i] = _mm512_unpacklo_epi8(reinterpret_cast<__m512i const &>(matrix[2 * i]),
                                       reinterpret_cast<__m512i const &>(matrix[2 * i + 1]));
        tmp1[i + 32] = _mm512_unpackhi_epi8(reinterpret_cast<__m512i const &>(matrix[2 * i]),
                                            reinterpret_cast<__m512i const &>(matrix[2 * i + 1]));
    }

    // Step 2: Unpack 16-bit operands.
    __m512i tmp2[64];
    for (int i = 0; i < 32; ++i)
    {
        tmp2[i] = _mm512_unpacklo_epi16(tmp1[2 * i], tmp1[2 * i + 1]);
        tmp2[i + 32] = _mm512_unpackhi_epi16(tmp1[2 * i], tmp1[2 * i + 1]);
    }
    // Step 3: Unpack 32-bit operands.
    for (int i = 0; i < 32; ++i)
    {
        tmp1[i] = _mm512_unpacklo_epi32(tmp2[2 * i], tmp2[2 * i + 1]);
        tmp1[i + 32] = _mm512_unpackhi_epi32(tmp2[2 * i], tmp2[2 * i + 1]);
    }
    // Step 4: Unpack 64-bit operands.
    for (int i = 0; i < 32; ++i)
    {
        tmp2[i] = _mm512_unpacklo_epi64(tmp1[2 * i], tmp1[2 * i + 1]);
        tmp2[i + 32] = _mm512_unpackhi_epi64(tmp1[2 * i], tmp1[2 * i + 1]);
    }

    // Step 5: Unpack 128-bit operands.
    // Helper function to emulate unpack of 128-bit across lanes using _mm512_permutex2var_epi64.
    auto _mm512_unpacklo_epi128 = [](__m512i const & a, __m512i const & b)
    {
        constexpr std::array<uint64_t, 8> lo_mask{0, 1, 8, 9, 2, 3, 10, 11};
        return _mm512_permutex2var_epi64(a, reinterpret_cast<__m512i const &>(*lo_mask.data()), b);
    };

    auto _mm521_unpackhi_epi128 = [](__m512i const & a, __m512i const & b)
    {
        constexpr std::array<uint64_t, 8> hi_mask{4, 5, 12, 13, 6, 7, 14, 15};
        return _mm512_permutex2var_epi64(a, reinterpret_cast<__m512i const &>(*hi_mask.data()), b);
    };

    for (int i = 0; i < 32; ++i)
    {
        tmp1[i] = _mm512_unpacklo_epi128(tmp2[2 * i], tmp2[2 * i + 1]);
        tmp1[i + 32] = _mm521_unpackhi_epi128(tmp2[2 * i], tmp2[2 * i + 1]);
    }
    // Step 6: Unpack 128-bit operands.
    // Helper function to emulate unpack of 256-bit across lanes using _mm512_shuffle_i64x2.
    auto _mm512_unpacklo_epi256 = [](__m512i const & a, __m512i const & b)
    {
        return _mm512_shuffle_i64x2(a, b, 0b0100'0100);
    };

    auto _mm521_unpackhi_epi256 = [](__m512i const & a, __m512i const & b)
    {
        return _mm512_shuffle_i64x2(a, b, 0b1110'1110);
    };

    // A look-up table to place the final transposed vector to the correct position in the
    // original matrix.
    // clang-format off
    constexpr std::array<uint32_t, 64> reverse_idx_mask{
    // 00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15
        0, 16,  8, 24,  4, 20, 12, 28,  2, 18, 10, 26,  6, 22, 14, 30,
    // 16  17  18  19  20  21  22  23  24  25  26  27  28  29  30  31
        1, 17,  9, 25,  5, 21, 13, 29,  3, 19, 11, 27,  7, 23, 15, 31,
    // 32  33  34  35  36  37  38  39  40  41  42  43  44  45  46  47
       32, 48, 40, 56, 36, 52, 44, 60, 34, 50, 42, 58, 38, 54, 46, 62,
    // 48  49  50  51  52  53  54  55  56  57  58  59  60  61  62  63
       33, 49, 41, 57, 37, 53, 45, 61, 35, 51, 43, 59, 39, 55, 47, 63};
    // clang-format on

    for (int i = 0; i < 32; ++i)
    {
        int const idx = i * 2;
        matrix[reverse_idx_mask[idx]] = reinterpret_cast<simd_t>(_mm512_unpacklo_epi256(tmp1[idx], tmp1[idx + 1]));
        matrix[reverse_idx_mask[idx + 1]] = reinterpret_cast<simd_t>(_mm521_unpackhi_epi256(tmp1[idx], tmp1[idx + 1]));
    }
}
#    endif // defined(__AVX512BW__)

template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed_avx512(source_simd_t const & src)
{
    __m512i const & tmp = reinterpret_cast<__m512i const &>(src);
    if constexpr (simd_traits<source_simd_t>::length == 64) // cast from epi8 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 32) // to epi16
            return reinterpret_cast<target_simd_t>(_mm512_cvtepi8_epi16(_mm512_castsi512_si256(tmp)));
        if constexpr (simd_traits<target_simd_t>::length == 16) // to epi32
            return reinterpret_cast<target_simd_t>(_mm512_cvtepi8_epi32(_mm512_castsi512_si128(tmp)));
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi64
            return reinterpret_cast<target_simd_t>(_mm512_cvtepi8_epi64(_mm512_castsi512_si128(tmp)));
    }
    else if constexpr (simd_traits<source_simd_t>::length == 32) // cast from epi16 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 16) // to epi32
            return reinterpret_cast<target_simd_t>(_mm512_cvtepi16_epi32(_mm512_castsi512_si256(tmp)));
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi64
            return reinterpret_cast<target_simd_t>(_mm512_cvtepi16_epi64(_mm512_castsi512_si128(tmp)));
    }
    else // cast from epi32 to epi64
    {
        static_assert(simd_traits<source_simd_t>::length == 16, "Expected 32 bit scalar type.");
        return reinterpret_cast<target_simd_t>(_mm512_cvtepi32_epi64(_mm512_castsi512_si256(tmp)));
    }
}

template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned_avx512(source_simd_t const & src)
{
    __m512i const & tmp = reinterpret_cast<__m512i const &>(src);
    if constexpr (simd_traits<source_simd_t>::length == 64) // cast from epi8 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 32) // to epi16
            return reinterpret_cast<target_simd_t>(_mm512_cvtepu8_epi16(_mm512_castsi512_si256(tmp)));
        if constexpr (simd_traits<target_simd_t>::length == 16) // to epi32
            return reinterpret_cast<target_simd_t>(_mm512_cvtepu8_epi32(_mm512_castsi512_si128(tmp)));
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi64
            return reinterpret_cast<target_simd_t>(_mm512_cvtepu8_epi64(_mm512_castsi512_si128(tmp)));
    }
    else if constexpr (simd_traits<source_simd_t>::length == 32) // cast from epi16 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 16) // to epi32
            return reinterpret_cast<target_simd_t>(_mm512_cvtepu16_epi32(_mm512_castsi512_si256(tmp)));
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi64
            return reinterpret_cast<target_simd_t>(_mm512_cvtepu16_epi64(_mm512_castsi512_si128(tmp)));
    }
    else // cast from epi32 to epi64
    {
        static_assert(simd_traits<source_simd_t>::length == 16, "Expected 32 bit scalar type.");
        return reinterpret_cast<target_simd_t>(_mm512_cvtepu32_epi64(_mm512_castsi512_si256(tmp)));
    }
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_half_avx512(simd_t const & src)
{
    return reinterpret_cast<simd_t>(
        _mm512_castsi256_si512(_mm512_extracti64x4_epi64(reinterpret_cast<__m512i const &>(src), index)));
}

#    if defined(__AVX512DQ__)
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_avx512(simd_t const & src)
{
    return reinterpret_cast<simd_t>(
        _mm512_castsi128_si512(_mm512_extracti64x2_epi64(reinterpret_cast<__m512i const &>(src), index)));
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_avx512(simd_t const & src)
{
    __m512i tmp = reinterpret_cast<__m512i const &>(src);

    // for uneven index exchange higher 64 bits with lower 64 bits for each 128 bit lane.
    if constexpr (index % 2 == 1)
        tmp = _mm512_shuffle_epi32(tmp, static_cast<_MM_PERM_ENUM>(0b0100'1110)); // := [1, 0, 3, 2].

    return reinterpret_cast<simd_t>(_mm512_castsi128_si512(_mm512_extracti64x2_epi64(tmp, index / 2)));
}
#    endif // defined(__AVX512DQ__)

} // namespace seqan3::detail

#endif // __AVX512F__
