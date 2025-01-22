// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides algorithms to modify seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cassert>
#include <concepts>
#include <utility>

#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/detail/builtin_simd.hpp>
#include <seqan3/utility/simd/detail/simd_algorithm_avx2.hpp>
#include <seqan3/utility/simd/detail/simd_algorithm_avx512.hpp>
#include <seqan3/utility/simd/detail/simd_algorithm_sse4.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

namespace seqan3::detail
{

//!\brief Helper function for seqan3::simd::fill.
//!\ingroup utility_simd
template <simd::simd_concept simd_t, size_t... I>
constexpr simd_t fill_impl(typename simd_traits<simd_t>::scalar_type const scalar, std::index_sequence<I...>) noexcept
{
    return simd_t{((void)I, scalar)...};
}

//!\brief Helper function for seqan3::simd::iota.
//!\ingroup utility_simd
template <simd::simd_concept simd_t, typename scalar_t, scalar_t... I>
constexpr simd_t iota_impl(scalar_t const offset, std::integer_sequence<scalar_t, I...>)
{
    return simd_t{static_cast<scalar_t>(offset + I)...};
}

/*!\brief Helper function to extract a part of the given simd vector.
 * \ingroup utility_simd
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
 * \ingroup utility_simd
 * \tparam target_simd_t The target simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \tparam source_simd_t The source simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \param[in] src The source to upcast into `target_simd_t`.
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
 * \ingroup utility_simd
 * \tparam target_simd_t The target simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \tparam source_simd_t The source simd type; must model seqan3::simd::simd_concept and must be a native builtin simd
 *                       type.
 * \param[in] src The source to upcast into `target_simd_t`.
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

/*!\brief Extracts one half of the given simd vector and stores the result in the lower half of the target vector.
 * \ingroup utility_simd
 * \tparam index An index value in the range of [0, 1].
 * \tparam simd_t The simd type.
 * \param src The source to extract the half from.
 * \returns A simd vector with the lower bits set with the respective half from `src`
 *          If the simd vector contains less than 2 elements, the unchanged source will be returned.
 *
 * \details
 *
 * Only the first simd length / 2 elements are defined.
 * The value of the remaining elements is implementation defined.
 *
 * \include test/snippet/utility/simd/simd_extract.cpp
 *
 * Example operation for SSE4:
 *
 * ```
 * i := index * 64
 * dst[63:0] := src[i+63:i]
 * ```
 */
template <uint8_t index, simd::simd_concept simd_t>
constexpr simd_t extract_half(simd_t const & src)
{
    static_assert(index < 2, "The index must be in the range of [0, 1]");

    return detail::extract_impl<2>(src, index);
}

//!\cond
template <uint8_t index, simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> && detail::is_native_builtin_simd_v<simd_t>
constexpr simd_t extract_half(simd_t const & src)
{
    static_assert(index < 2, "The index must be in the range of [0, 1]");

    if constexpr (simd_traits<simd_t>::length < 2) // In case there are less elements available return unchanged value.
        return src;
    else if constexpr (simd_traits<simd_t>::max_length == 16) // SSE4
        return detail::extract_half_sse4<index>(src);
    else if constexpr (simd_traits<simd_t>::max_length == 32) // AVX2
        return detail::extract_half_avx2<index>(src);
    else if constexpr (simd_traits<simd_t>::max_length == 64) // AVX512
        return detail::extract_half_avx512<index>(src);
    else // Anything else
        return detail::extract_impl<2>(src, index);
}
//!\endcond

/*!\brief Extracts one quarter of the given simd vector and stores it in the lower quarter of the target vector.
 * \ingroup utility_simd
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
 * \include test/snippet/utility/simd/simd_extract.cpp
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
    requires detail::is_builtin_simd_v<simd_t> && detail::is_native_builtin_simd_v<simd_t>
constexpr simd_t extract_quarter(simd_t const & src)
{
    static_assert(index < 4, "The index must be in the range of [0, 1, 2, 3]");

    if constexpr (simd_traits<simd_t>::length < 4) // In case there are less elements available return unchanged value.
        return src;
    else if constexpr (simd_traits<simd_t>::max_length == 16) // SSE4
        return detail::extract_quarter_sse4<index>(src);
    else if constexpr (simd_traits<simd_t>::max_length == 32) // AVX2
        return detail::extract_quarter_avx2<index>(src);
#if defined(__AVX512DQ__)
    else if constexpr (simd_traits<simd_t>::max_length == 64) // AVX512
        return detail::extract_quarter_avx512<index>(src);
#endif   // defined(__AVX512DQ__)
    else // Anything else
        return detail::extract_impl<4>(src, index);
}
//!\endcond

/*!\brief Extracts one eighth of the given simd vector and stores it in the lower eighth of the target vector.
 * \ingroup utility_simd
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
 * \include test/snippet/utility/simd/simd_extract.cpp
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
    requires detail::is_builtin_simd_v<simd_t> && detail::is_native_builtin_simd_v<simd_t>
constexpr simd_t extract_eighth(simd_t const & src)
{
    static_assert(index < 8, "The index must be in the range of [0, 1, 2, 3, 4, 5, 6, 7]");

    if constexpr (simd_traits<simd_t>::length < 8) // In case there are less elements available return unchanged value.
        return src;
    else if constexpr (simd_traits<simd_t>::max_length == 16) // SSE4
        return detail::extract_eighth_sse4<index>(src);
    else if constexpr (simd_traits<simd_t>::max_length == 32) // AVX2
        return detail::extract_eighth_avx2<index>(src);
#if defined(__AVX512DQ__)
    else if constexpr (simd_traits<simd_t>::max_length == 64) // AVX512
        return detail::extract_eighth_avx512<index>(src);
#endif   // defined(__AVX512DQ__)
    else // Anything else
        return detail::extract_impl<8>(src, index);
}
//!\endcond

//!\cond
template <simd::simd_concept simd_t>
constexpr void transpose(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    std::array<simd_t, simd_traits<simd_t>::length> tmp{};

    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            tmp[j][i] = matrix[i][j];

    std::swap(tmp, matrix);
}
//!\endcond
} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\brief Fills a seqan3::simd::simd_type vector with a scalar value.
 * \ingroup utility_simd
 * \tparam    simd_t The simd type which satisfies seqan3::simd::simd_concept.
 * \param[in] scalar The scalar value to fill the seqan3::simd::simd_type vector.
 *
 * \details
 *
 * \include test/snippet/utility/simd/fill.cpp
 */
template <simd::simd_concept simd_t>
constexpr simd_t fill(typename simd_traits<simd_t>::scalar_type const scalar) noexcept
{
    constexpr size_t length = simd_traits<simd_t>::length;
    return detail::fill_impl<simd_t>(scalar, std::make_index_sequence<length>{});
}

/*!\brief Fills a seqan3::simd::simd_type vector with the scalar values offset, offset+1, offset+2, ...
 * \ingroup utility_simd
 * \tparam    simd_t The simd type which satisfies seqan3::simd::simd_concept.
 * \param[in] offset The scalar offset to fill the seqan3::simd::simd_type vector
 *
 * \details
 *
 * \include test/snippet/utility/simd/iota.cpp
 */
template <simd::simd_concept simd_t>
constexpr simd_t iota(typename simd_traits<simd_t>::scalar_type const offset)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    using scalar_type = typename simd_traits<simd_t>::scalar_type;
    return detail::iota_impl<simd_t>(offset, std::make_integer_sequence<scalar_type, length>{});
}

/*!\brief Load simd_t size bits of integral data from memory.
 * \ingroup utility_simd
 * \tparam    simd_t   The simd type; must model seqan3::simd::simd_concept.
 * \param[in] mem_addr The memory address to load from. Does not need to be aligned on any particular boundary.
 *
 * \details
 *
 * \include test/snippet/utility/simd/simd_load.cpp
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
    requires detail::is_builtin_simd_v<simd_t> && detail::is_native_builtin_simd_v<simd_t>
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

/*!\brief Store simd_t size bits of integral data into memory.
 * \ingroup utility_simd
 * \tparam    simd_t   The simd type; must model seqan3::simd::simd_concept.
 * \param[in] mem_addr The memory address to write to. Does not need to be aligned on any particular boundary.
 * \param[in] simd_vec The simd vector to read the data from.
 *
 * \details
 *
 * \include test/snippet/utility/simd/simd_store.cpp
 */
template <simd::simd_concept simd_t>
constexpr void store(void * mem_addr, simd_t const & simd_vec)
{
    assert(mem_addr != nullptr);
    using scalar_t = typename simd_traits<simd_t>::scalar_type;

    for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        *(static_cast<scalar_t *>(mem_addr) + i) = simd_vec[i];
}

//!\cond
template <simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> && detail::is_native_builtin_simd_v<simd_t>
constexpr void store(void * mem_addr, simd_t const & simd_vec)
{
    assert(mem_addr != nullptr);

    if constexpr (simd_traits<simd_t>::max_length == 16)
        detail::store_sse4<simd_t>(mem_addr, simd_vec);
    else if constexpr (simd_traits<simd_t>::max_length == 32)
        detail::store_avx2<simd_t>(mem_addr, simd_vec);
    else if constexpr (simd_traits<simd_t>::max_length == 64)
        detail::store_avx512<simd_t>(mem_addr, simd_vec);
    else
        static_assert(simd_traits<simd_t>::max_length >= 16 && simd_traits<simd_t>::max_length <= 64,
                      "Unsupported simd type.");
}
//!\endcond

/*!\brief Transposes the given simd vector matrix.
 * \ingroup utility_simd
 * \tparam simd_t The simd vector type; must model seqan3::simd::simd_concept and must be a simd built-in type.
 * \param[in,out] matrix The matrix that is transposed in place.
 *
 * \details
 *
 * \include test/snippet/utility/simd/simd_transpose.cpp
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
    detail::transpose(matrix);
}

//!\cond
// Implementation for seqan builtin simd.
template <simd::simd_concept simd_t>
    requires detail::is_builtin_simd_v<simd_t> && detail::is_native_builtin_simd_v<simd_t>
          && (simd_traits<simd_t>::max_length == simd_traits<simd_t>::length)
constexpr void transpose(std::array<simd_t, simd_traits<simd_t>::length> & matrix)
{
    if constexpr (simd_traits<simd_t>::length == 16) // SSE4 implementation
        detail::transpose_matrix_sse4(matrix);
    else if constexpr (simd_traits<simd_t>::length == 32) // AVX2 implementation
        detail::transpose_matrix_avx2(matrix);
#if defined(__AVX512BW__)                                 // Requires byte-word extension of AVX512 instruction set.
    else if constexpr (simd_traits<simd_t>::length == 64) // AVX512 implementation
        detail::transpose_matrix_avx512(matrix);
#endif // defined(__AVX512BW__)
    else
        detail::transpose(matrix);
}
//!\endcond

/*!\brief Upcasts the given vector into the target vector using sign extension of packed values.
 * \ingroup utility_simd
 * \tparam simd_t The simd type; must model seqan3::simd::simd_concept and must be a builtin simd type.
 *
 * \details
 *
 * \include test/snippet/utility/simd/simd_upcast.cpp
 */
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
constexpr target_simd_t upcast(source_simd_t const & src)
{
    static_assert(
        simd_traits<target_simd_t>::length <= simd_traits<source_simd_t>::length,
        "The length of the target simd type must be greater or equal than the length of the source simd type.");

    target_simd_t tmp{};
    for (unsigned i = 0; i < simd_traits<target_simd_t>::length; ++i)
        tmp[i] = static_cast<typename simd_traits<target_simd_t>::scalar_type>(src[i]);

    return tmp;
}

//!\cond
template <simd::simd_concept target_simd_t, simd::simd_concept source_simd_t>
    requires detail::is_builtin_simd_v<target_simd_t> && detail::is_builtin_simd_v<source_simd_t>
          && detail::is_native_builtin_simd_v<source_simd_t>
constexpr target_simd_t upcast(source_simd_t const & src)
{
    static_assert(
        simd_traits<target_simd_t>::length <= simd_traits<source_simd_t>::length,
        "The length of the target simd type must be greater or equal than the length of the source simd type.");

    if constexpr (simd_traits<source_simd_t>::length == simd_traits<target_simd_t>::length)
    {
        static_assert(simd_traits<target_simd_t>::max_length == simd_traits<source_simd_t>::max_length,
                      "Target vector has a different byte size.");
        return reinterpret_cast<target_simd_t>(src); // Same packing so we do not cast.
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

} // namespace simd

} // namespace seqan3
