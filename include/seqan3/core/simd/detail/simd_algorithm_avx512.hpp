// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides specific algorithm implementations for AVX512 instruction set.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/builtin_simd_intrinsics.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

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

/*!\copydoc seqan3::detail::extract_halve
 * \attention This is the implementation for AVX512 intrinsics.
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_halve_avx512(simd_t const & src);

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

}

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

// TODO: not implemented and used yet, if you implement it don't forget to add it to seqan3::simd::transpose
template <simd::simd_concept simd_t>
inline void transpose_matrix_avx512(std::array<simd_t, simd_traits<simd_t>::length> & matrix);

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

// TODO: not implemented and used yet, if you implement it don't forget to add it to seqan3::detail::extract_halve
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_halve_avx512(simd_t const & src);

// TODO: not implemented and used yet, if you implement it don't forget to add it to seqan3::detail::extract_quarter
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter_avx512(simd_t const & src);

// TODO: not implemented and used yet, if you implement it don't forget to add it to seqan3::detail::extract_eighth
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth_avx512(simd_t const & src);

} // namespace seqan3::detail

#endif // __AVX512F__
