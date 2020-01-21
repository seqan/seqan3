// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides specific algorithm implementations for SSE4 instruction set.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/builtin_simd_intrinsics.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

//-----------------------------------------------------------------------------
// forward declare sse4 simd algorithms that use sse4 intrinsics
//-----------------------------------------------------------------------------

namespace seqan3::detail
{
/*!\copydoc seqan3::simd::load
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <simd::simd_concept simd_t>
constexpr simd_t load_sse4(void const * mem_addr);

/*!\copydoc seqan3::simd::transpose
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <simd::simd_concept simd_t>
inline void transpose_matrix_sse4(std::array<simd_t, simd_traits<simd_t>::length> & matrix);

/*!\copydoc seqan3::detail::upcast_signed
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed_sse4(source_simd_t const & src);

/*!\copydoc seqan3::detail::upcast_unsigned
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned_sse4(source_simd_t const & src);

/*!\copydoc seqan3::detail::extract_halve
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_halve_sse4(simd_t const & src);

/*!\copydoc seqan3::detail::extract_quarter
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_sse4(simd_t const & src);

/*!\copydoc seqan3::detail::extract_eighth
 * \attention This is the implementation for SSE4 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_sse4(simd_t const & src);

}

//-----------------------------------------------------------------------------
// implementation
//-----------------------------------------------------------------------------

#ifdef __SSE4_2__

namespace seqan3::detail
{

template <simd::simd_concept simd_t>
constexpr simd_t load_sse4(void const * mem_addr)
{
    return reinterpret_cast<simd_t>(_mm_loadu_si128(reinterpret_cast<__m128i const *>(mem_addr)));
}

template <simd::simd_concept simd_t>
inline void transpose_matrix_sse4(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    static_assert(simd_traits<simd_t>::length == simd_traits<simd_t>::max_length, "Expects byte scalar type.");
    static_assert(is_native_builtin_simd_v<simd_t>, "The passed simd vector is not a native SSE4 simd vector type.");
    static_assert(is_builtin_simd_v<simd_t>, "The passed simd vector is not a builtin vector type.");

    // we need a look-up table to reverse the lowest 4 bits
    // in order to place the permute the transposed rows
    constexpr std::array<char, 16> bit_reverse{0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15};

    // transpose a 16x16 byte matrix
    //
    // matrix =
    // A0 A1 A2 ... Ae Af
    // B0 B1 B2 ... Be Bf
    // ...
    // P0 P1 P2 ... Pe Pf
    __m128i tmp1[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi8(reinterpret_cast<__m128i &>(matrix[2*i]),
                                      reinterpret_cast<__m128i &>(matrix[2*i+1]));
        tmp1[i+8] = _mm_unpackhi_epi8(reinterpret_cast<__m128i &>(matrix[2*i]),
                                      reinterpret_cast<__m128i &>(matrix[2*i+1]));
    }
    // tmp1[0]  = A0 B0 A1 B1 ... A7 B7
    // tmp1[1]  = C0 D0 C1 D1 ... C7 D7
    // ...
    // tmp1[7]  = O0 P0 O1 P1 ... O7 P7
    // tmp1[8]  = A8 B8 A9 B9 ... Af Bf
    // ...
    // tmp1[15] = O8 P8 O9 P9 ... Of Pf
    __m128i tmp2[16];
    for (int i = 0; i < 8; ++i)
    {
        tmp2[i]   = _mm_unpacklo_epi16(tmp1[2*i], tmp1[2*i+1]);
        tmp2[i+8] = _mm_unpackhi_epi16(tmp1[2*i], tmp1[2*i+1]);
    }
    // tmp2[0]  = A0 B0 C0 D0 ... A3 B3 C3 D3
    // tmp2[1]  = E0 F0 G0 H0 ... E3 F3 G3 H3
    // ...
    // tmp2[3]  = M0 N0 O0 P0 ... M3 N3 O3 P3
    // tmp2[4]  = A8 B8 C8 D8 ... Ab Bb Cb Db
    // ...
    // tmp2[7]  = M8 N8 O8 P8 ... Mb Nb Ob Pb
    // tmp2[8]  = A4 B4 C4 D4 ... A7 B7 C7 D7
    // ..
    // tmp2[12] = Ac Bc Cc Dc ... Af Bf Cf Df
    // ...
    // tmp2[15] = Mc Nc Oc Pc ... Mf Nf Of Pf
    for (int i = 0; i < 8; ++i)
    {
        tmp1[i]   = _mm_unpacklo_epi32(tmp2[2*i], tmp2[2*i+1]);
        tmp1[i+8] = _mm_unpackhi_epi32(tmp2[2*i], tmp2[2*i+1]);
    }
    // tmp1[0]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    // ...
    // tmp1[4]  = A0 B0 .... H0 A1 B1 .... H1
    // tmp1[1]  = I0 J0 .... P0 I1 J1 .... P1
    for (int i = 0; i < 8; ++i)
    {
        matrix[bit_reverse[i]]   = reinterpret_cast<simd_t>(_mm_unpacklo_epi64(tmp1[2*i], tmp1[2*i+1]));
        matrix[bit_reverse[i+8]] = reinterpret_cast<simd_t>(_mm_unpackhi_epi64(tmp1[2*i], tmp1[2*i+1]));
    }
}

template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed_sse4(source_simd_t const & src)
{
    if constexpr (simd_traits<source_simd_t>::length == 16) // cast from epi8 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi16
            return reinterpret_cast<target_simd_t>(_mm_cvtepi8_epi16(reinterpret_cast<__m128i const &>(src)));
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi32
            return reinterpret_cast<target_simd_t>(_mm_cvtepi8_epi32(reinterpret_cast<__m128i const &>(src)));
        if constexpr (simd_traits<target_simd_t>::length == 2) // to epi64
            return reinterpret_cast<target_simd_t>(_mm_cvtepi8_epi64(reinterpret_cast<__m128i const &>(src)));
    }
    else if constexpr (simd_traits<source_simd_t>::length == 8) // cast from epi16 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi32
            return reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi32(reinterpret_cast<__m128i const &>(src)));
        if constexpr (simd_traits<target_simd_t>::length == 2) // to epi64
            return reinterpret_cast<target_simd_t>(_mm_cvtepi16_epi64(reinterpret_cast<__m128i const &>(src)));
    }
    else  // cast from epi32 to epi64
    {
        static_assert(simd_traits<source_simd_t>::length == 4, "Expected 32 bit scalar type.");
        return reinterpret_cast<target_simd_t>(_mm_cvtepi32_epi64(reinterpret_cast<__m128i const &>(src)));
    }
}

template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned_sse4(source_simd_t const & src)
{
    if constexpr (simd_traits<source_simd_t>::length == 16) // cast from epi8 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 8) // to epi16
            return reinterpret_cast<target_simd_t>(_mm_cvtepu8_epi16(reinterpret_cast<__m128i const &>(src)));
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi32
            return reinterpret_cast<target_simd_t>(_mm_cvtepu8_epi32(reinterpret_cast<__m128i const &>(src)));
        if constexpr (simd_traits<target_simd_t>::length == 2) // to epi64
            return reinterpret_cast<target_simd_t>(_mm_cvtepu8_epi64(reinterpret_cast<__m128i const &>(src)));
    }
    else if constexpr (simd_traits<source_simd_t>::length == 8) // cast from epi16 ...
    {
        if constexpr (simd_traits<target_simd_t>::length == 4) // to epi32
            return reinterpret_cast<target_simd_t>(_mm_cvtepu16_epi32(reinterpret_cast<__m128i const &>(src)));
        if constexpr (simd_traits<target_simd_t>::length == 2) // to epi64
            return reinterpret_cast<target_simd_t>(_mm_cvtepu16_epi64(reinterpret_cast<__m128i const &>(src)));
    }
    else  // cast from epi32 to epi64
    {
        static_assert(simd_traits<source_simd_t>::length == 4, "Expected 32 bit scalar type.");
        return reinterpret_cast<target_simd_t>(_mm_cvtepu32_epi64(reinterpret_cast<__m128i const &>(src)));
    }
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_halve_sse4(simd_t const & src)
{
    return reinterpret_cast<simd_t>(_mm_srli_si128(reinterpret_cast<__m128i const &>(src), (index) << 3));
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_sse4(simd_t const & src)
{
    return reinterpret_cast<simd_t>(_mm_srli_si128(reinterpret_cast<__m128i const &>(src), index << 2));
}

template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_sse4(simd_t const & src)
{
    return reinterpret_cast<simd_t>(_mm_srli_si128(reinterpret_cast<__m128i const &>(src), index << 1));
}

} // namespace seqan3::detail

#endif // __SSE4_2__
