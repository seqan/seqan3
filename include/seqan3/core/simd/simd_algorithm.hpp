// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithms to modify seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <utility>

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/builtin_simd.hpp>
#include <seqan3/core/simd/detail/simd_algorithm_sse4.hpp>
#include <seqan3/core/simd/detail/simd_algorithm_avx2.hpp>
#include <seqan3/core/simd/detail/simd_algorithm_avx512.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

//!\brief Helper function for seqan3::simd::fill.
//!\ingroup simd
template <simd::simd_concept simd_t, size_t... I>
constexpr simd_t fill_impl(typename simd_traits<simd_t>::scalar_type const scalar, std::index_sequence<I...>) noexcept
{
    return simd_t{((void)I, scalar)...};
}

//!\brief Helper function for seqan3::simd::iota.
//!\ingroup simd
template <simd::simd_concept simd_t, typename scalar_t, scalar_t... I>
constexpr simd_t iota_impl(scalar_t const offset, std::integer_sequence<scalar_t, I...>)
{
    return simd_t{static_cast<scalar_t>(offset + I)...};
}

/*!\brief Helper function to extract a part of the given simd vector.
 * \ingroup simd
 * \tparam divisor The divisor to select the chunk size.
 * \tparam simd_t  The simd type; must model seqan3::simd::simd_concept.
 *
 * \param[in] src  The source vector to extract from.
 * \param[in] mask The control mask to select which chunk is extracted.
 * \returns The destination vector containing the extracted part.
 *
 * \details
 *
 * Extracts the specified part of the source simd vector and stores it in the first chunk starting at offset 0 in
 * the destination vector.
 */
template <size_t divisor, simd_concept simd_t>
constexpr simd_t extract_impl(simd_t const & src, uint8_t const mask)
{
    simd_t dst{};
    constexpr size_t chunk = simd_traits<simd_t>::length / divisor;
    size_t offset = chunk * mask;
    for (size_t i = 0; i < chunk; ++i)
        dst[i] = src[i + offset];

    return dst;
}

/*!\brief Upcasts the given vector into the target vector using signed extension of packed values.
 * \tparam target_simd_t The target simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \tparam source_simd_t The source simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \param[in] src The source to upcast into `target_simd_t`.
 * \ingroup simd
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_signed(source_simd_t const & src)
{
    static_assert(simd_traits<target_simd_t>::max_length == simd_traits<source_simd_t>::max_length,
                  "Target vector has different byte size.");

    if constexpr (simd_traits<source_simd_t>::max_length == 16) // SSE4
        return upcast_signed_sse4<target_simd_t>(src);
    else if constexpr (simd_traits<source_simd_t>::max_length == 32) // AVX2
        return upcast_signed_avx2<target_simd_t>(src);
    else if constexpr (simd_traits<source_simd_t>::max_length == 64) // AVX512
        return upcast_signed_avx512<target_simd_t>(src);
    else
        static_assert(simd_traits<source_simd_t>::max_length <= 32, "simd type is not supported.");
}

/*!\brief Upcasts the given vector into the target vector using unsigned extension of packed values.
 * \tparam target_simd_t The target simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \tparam source_simd_t The source simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \param[in] src The source to upcast into `target_simd_t`.
 * \ingroup simd
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast_unsigned(source_simd_t const & src)
{
    static_assert(simd_traits<target_simd_t>::max_length == simd_traits<source_simd_t>::max_length,
                  "Target vector has different byte size.");

    if constexpr (simd_traits<source_simd_t>::max_length == 16) // SSE4
        return upcast_unsigned_sse4<target_simd_t>(src);
    else if constexpr (simd_traits<source_simd_t>::max_length == 32) // AVX2
        return upcast_unsigned_avx2<target_simd_t>(src);
    else if constexpr (simd_traits<source_simd_t>::max_length == 64) // AVX512
        return upcast_unsigned_avx512<target_simd_t>(src);
    else
        static_assert(simd_traits<source_simd_t>::max_length <= 32, "simd type is not supported.");
}

/*!\brief Extracts one halve of the given simd vector and stores the result in the lower halve of the target vector.
 * \ingroup simd
 * \tparam index An index value in the range of [0, 1].
 * \tparam simd_t The simd type.
 * \param src The source to extract the halve from.
 * \returns A simd vector with the lower bits set with the respective halve from `src`
 *          If the simd vector contains less than 2 elements, the unchanged source will be returned.
 *
 * \details
 *
 * Only the first simd length / 2 elements are defined.
 * The value of the remaining elements is implementation defined.
 *
 * \include test/snippet/core/simd/simd_extract.cpp
 *
 * Example operation for SSE4:
 *
 * ```
 * i := index * 64
 * dst[63:0] := src[i+63:i]
 * ```
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_halve(simd_t const & src)
{
    static_assert(index < 2, "The index must be in the range of [0, 1]");

    return detail::extract_impl<2>(src, index);
}

//!\cond
template <uint8_t index, simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> &&
             detail::is_native_builtin_simd_v<simd_t> &&
             simd_traits<simd_t>::max_length == 16
constexpr simd_t extract_halve(simd_t const & src)
{
    static_assert(index < 2, "The index must be in the range of [0, 1]");

    if constexpr (simd_traits<simd_t>::length < 2) // In case there are less elements available return unchanged value.
        return src;
    else // if constexpr (simd_traits<simd_t>::max_length == 16) // SSE4
        return detail::extract_halve_sse4<index>(src);
}
//!\endcond

/*!\brief Extracts one quarter of the given simd vector and stores it in the lower quarter of the target vector.
 * \ingroup simd
 * \tparam index An index value in the range of [0, 1, 2, 3].
 * \tparam simd_t The simd type.
 * \param src The source to extract the quarter from.
 * \returns A simd vector with the lower bits set with the respective quarter from `src`.
 *          If the simd vector contains less than 4 elements, the unchanged source will be returned.
 *
 * \details
 *
 * Only the first simd length / 4 elements are defined.
 * The value of the remaining elements is implementation defined.
 *
 * \include test/snippet/core/simd/simd_extract.cpp
 *
 * Example operation for SSE4:
 *
 * ```
 * i := index * 32
 * dst[31:0] := src[i+31:i]
 * ```
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_quarter(simd_t const & src)
{
    static_assert(index < 4, "The index must be in the range of [0, 1, 2, 3]");

    return detail::extract_impl<4>(src, index);
}

//!\cond
template <uint8_t index, simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> &&
             detail::is_native_builtin_simd_v<simd_t> &&
             simd_traits<simd_t>::max_length == 16
constexpr simd_t extract_quarter(simd_t const & src)
{
    static_assert(index < 4, "The index must be in the range of [0, 1, 2, 3]");

    if constexpr (simd_traits<simd_t>::length < 4) // In case there are less elements available return unchanged value.
        return src;
    else // if constexpr (simd_traits<simd_t>::max_length == 16) // SSE4
        return detail::extract_quarter_sse4<index>(src);
}
//!\endcond

/*!\brief Extracts one eighth of the given simd vector and stores it in the lower eighth of the target vector.
 * \ingroup simd
 * \tparam index An index value in the range of [0, 1, 2, 3, 4, 5, 6, 7].
 * \tparam simd_t The simd type.
 * \param src The source to extract the eighth from.
 * \returns A simd vector with the lower bits set with the respective eighth from `src`
 *          If the simd vector contains less than 8 elements, the unchanged source will be returned.
 *
 * \details
 *
 * Only the first simd length / 8 elements are defined.
 * The value of the remaining elements is implementation defined.
 *
 * \include test/snippet/core/simd/simd_extract.cpp
 *
 * Example operation for SSE4:
 *
 * ```
 * i := index * 16
 * dst[15:0] := src[i+15:i]
 * ```
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_eighth(simd_t const & src)
{
    return detail::extract_impl<8>(src, index);
}

//!\cond
template <uint8_t index, simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> &&
             detail::is_native_builtin_simd_v<simd_t> &&
             simd_traits<simd_t>::max_length == 16
constexpr simd_t extract_eighth(simd_t const & src)
{
    static_assert(index < 8, "The index must be in the range of [0, 1, 2, 3, 4, 5, 6, 7]");

    if constexpr (simd_traits<simd_t>::length < 8) // In case there are less elements available return unchanged value.
        return src;
    else // if constexpr (simd_traits<simd_t>::max_length == 16) // SSE4
        return detail::extract_eighth_sse4<index>(src);
}
//!\endcond

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\brief Fills a seqan3::simd::simd_type vector with a scalar value.
 * \tparam    simd_t The simd type which satisfies seqan3::simd::simd_concept.
 * \param[in] scalar The scalar value to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/fill.cpp
 */
template <simd::simd_concept simd_t>
constexpr simd_t fill(typename simd_traits<simd_t>::scalar_type const scalar) noexcept
{
    constexpr size_t length = simd_traits<simd_t>::length;
    return detail::fill_impl<simd_t>(scalar, std::make_index_sequence<length>{});
}

/*!\brief Fills a seqan3::simd::simd_type vector with the scalar values offset, offset+1, offset+2, ...
 * \tparam    simd_t The simd type which satisfies seqan3::simd::simd_concept.
 * \param[in] offset The scalar offset to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/iota.cpp
 */
template <simd::simd_concept simd_t>
constexpr simd_t iota(typename simd_traits<simd_t>::scalar_type const offset)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    using scalar_type = typename simd_traits<simd_t>::scalar_type;
    return detail::iota_impl<simd_t>(offset, std::make_integer_sequence<scalar_type, length>{});
}

/*!\brief Load simd_t size bits of integral data from memory.
 * \ingroup simd
 * \tparam    simd_t   The simd type; must model seqan3::simd::simd_concept.
 * \param[in] mem_addr The memory address to load from. Does not need to be aligned on any particular boundary.
 *
 * \details
 *
 * \include test/snippet/core/simd/simd_load.cpp
 */
template <simd::simd_concept simd_t>
constexpr simd_t load(void const * mem_addr)
{
    assert(mem_addr != nullptr);
    simd_t tmp{};

    for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        tmp[i] = *(static_cast<typename simd_traits<simd_t>::scalar_type const *>(mem_addr) + i);

    return tmp;
}

//!\cond
template <simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> &&
             detail::is_native_builtin_simd_v<simd_t>
constexpr simd_t load(void const * mem_addr)
{
    assert(mem_addr != nullptr);

    if constexpr (simd_traits<simd_t>::max_length == 16)
        return detail::load_sse4<simd_t>(mem_addr);
    else if constexpr (simd_traits<simd_t>::max_length == 32)
        return detail::load_avx2<simd_t>(mem_addr);
    else if constexpr (simd_traits<simd_t>::max_length == 64)
        return detail::load_avx512<simd_t>(mem_addr);
    else
        static_assert(simd_traits<simd_t>::max_length >= 16 && simd_traits<simd_t>::max_length <= 64,
                      "Unsupported simd type.");
}
//!\endcond

/*!\brief Transposes the given simd vector matrix.
 * \ingroup simd
 * \tparam simd_t The simd vector type; must model seqan3::simd::simd_concept and must be a simd built-in type.
 * \param[in,out] matrix The matrix that is transposed in place.
 *
 * \details
 *
 * \include test/snippet/core/simd/simd_transpose.cpp
 *
 * ### Exception
 *
 * Strong exception guarantee.
 *
 * ### Complexity
 *
 * Quadratic.
 */
template <simd::simd_concept simd_t>
constexpr void transpose(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    std::array<simd_t, simd_traits<simd_t>::length> tmp{};

    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            tmp[j][i] = matrix[i][j];

    std::swap(tmp, matrix);
}

//!\cond
template <simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> &&
             detail::is_native_builtin_simd_v<simd_t> &&
             simd_traits<simd_t>::max_length == 16 &&
             simd_traits<simd_t>::length == 16
constexpr void transpose(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    detail::transpose_matrix_sse4(matrix);
}
//!\endcond

/*!\brief Upcasts the given vector into the target vector using sign extension of packed values.
 * \ingroup simd
 * \tparam simd_t The simd type; must model seqan3::simd::simd_concept and must be a builtin simd type.
 *
 * \details
 *
 * \include test/snippet/core/simd/simd_upcast.cpp
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast(source_simd_t const & src)
{
    static_assert(simd_traits<target_simd_t>::length <= simd_traits<source_simd_t>::length,
                  "The length of the target simd type must be greater or equal than the length of the source simd type.");

    target_simd_t tmp{};
    for (unsigned i = 0; i < simd_traits<target_simd_t>::length; ++i)
        tmp[i] = static_cast<typename simd_traits<target_simd_t>::scalar_type>(src[i]);

    return tmp;
}

//!\cond
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
    requires detail::is_builtin_simd_v<target_simd_t> &&
             detail::is_builtin_simd_v<source_simd_t> &&
             detail::is_native_builtin_simd_v<source_simd_t>
constexpr target_simd_t upcast(source_simd_t const & src)
{
    static_assert(simd_traits<target_simd_t>::length <= simd_traits<source_simd_t>::length,
                  "The length of the target simd type must be greater or equal than the length of the source simd type.");

    if constexpr (simd_traits<source_simd_t>::length == simd_traits<target_simd_t>::length)
    {
        static_assert(simd_traits<target_simd_t>::max_length == simd_traits<source_simd_t>::max_length,
                    "Target vector has a different byte size.");
        return reinterpret_cast<target_simd_t>(src);  // Same packing so we do not cast.
    }
    else if constexpr (std::signed_integral<typename simd_traits<source_simd_t>::scalar_type>)
    {
        return detail::upcast_signed<target_simd_t>(src);
    }
    else
    {
        static_assert(std::unsigned_integral<typename simd_traits<source_simd_t>::scalar_type>,
                      "Expected unsigned scalar type.");
        return detail::upcast_unsigned<target_simd_t>(src);
    }
}
//!\endcond

} // inline namespace simd

} // namespace seqan3
