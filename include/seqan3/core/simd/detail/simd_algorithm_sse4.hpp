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

/*!\brief Implementation of seqan3::simd::unpack_hi for sse4.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[in] first   The vector whose values come before the `second`.
 * \param[in] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <Simd simd_t>
//!\cond
    requires simd_traits<std::remove_reference_t<simd_t>>::max_length == 16
//!\endcond
inline simd_t unpack_hi(simd_t const & first, simd_t const & second)
{
    constexpr size_t scalar_size = sizeof(typename simd_traits<simd_t>::scalar_type);

    if constexpr (scalar_size == 1)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi8(reinterpret_cast<__m128i const &>(first),
                                                          reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 2)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi16(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 4)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi32(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 8)
    {
        return reinterpret_cast<simd_t>(_mm_unpackhi_epi64(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }
}

/*!\brief Implementation of seqan3::simd::unpack_lo for sse4.
 * \tparam simd_t      The simd type; must model seqan3::simd::Simd.
 * \param[in] first   The vector whose values come before the `second`.
 * \param[in] second  The vector whose values come after the `first`.
 * \ingroup simd
 */
template <typename simd_t>
//!\cond
    requires simd_traits<simd_t>::max_length == 16
//!\endcond
inline simd_t unpack_lo(simd_t const & first, simd_t const & second)
{
    constexpr size_t scalar_size = sizeof(typename simd_traits<simd_t>::scalar_type);

    if constexpr (scalar_size == 1)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi8(reinterpret_cast<__m128i const &>(first),
                                                          reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 2)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi16(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 4)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi32(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else if constexpr (scalar_size == 8)
    {
        return reinterpret_cast<simd_t>(_mm_unpacklo_epi64(reinterpret_cast<__m128i const &>(first),
                                                           reinterpret_cast<__m128i const &>(second)));
    }
    else
    {
        static_assert(scalar_size <= 8 && !is_power_of_two(scalar_size),
                      "The targeted scalar size is not supported.");
    }
}

} // namespace seqan3::detail
