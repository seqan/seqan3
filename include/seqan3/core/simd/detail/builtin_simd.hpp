// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::builtin_simd, seqan3::detail::is_builtin_simd
 * and seqan3::simd::simd_traits<builtin_simd_t>.
 */

#pragma once

#include <type_traits>

#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/core/metafunction/transformation_trait_or.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/detail/builtin_simd_macros.hpp>
#include <seqan3/core/simd/detail/default_simd_max_length.hpp>

namespace seqan3::detail
{

/*!\brief seqan3::detail::builtin_simd is the class that holds the type of a
 * simd implementation called [vector extension]
 * (https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html)
 * (formerly known as "seqan simd" in seqan2).
 * \ingroup simd
 * \tparam scalar_type The underlying type of a simd vector
 * \tparam length_v    The number of packed values in a simd vector
 *
 * \include test/snippet/core/simd/detail/builtin_simd.cpp
 *
 * seqan3::detail::builtin_simd is basically defined as:
 *
 * ```
 * template <typename scalar_t, size_t length>
 * struct builtin_simd
 * {
 *    using type = scalar_t __attribute__ ((__vector_size__(length * sizeof(scalar_t))));
 * };
 * ```
 *
 * \attention This class itself only delegates to a [vector extension]
 * (https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html) type,
 * which is offered by the compiler as a builtin type.
 *
 * But gcc doesn't allow attributes on template arguments as this warning illustrates:
 * >   ignoring attributes applied to dependent type ‘scalar_t’ without an
 * >   associated declaration [-Wattributes]
 *
 * We use the following macro hell to define this explicitly for `[u]intX_t` types:
 * \include test/snippet/core/simd/detail/builtin_simd_macro.cpp
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html
 */
template <typename scalar_t, size_t length>
struct builtin_simd
#if SEQAN3_DOXYGEN_ONLY(1)0
{
    //!\brief The type of the builtin simd.
    using type = IMPLEMENTATION_DEFINED;
}
#endif
;

/*!\brief #SEQAN3_BUILTIN_SIMD_CLASS defines seqan3::detail::builtin_simd for
 * `[u]intX_t` types. Will be undefined immediately.
 * \ingroup simd
 * \param  scalar_t   same as seqan3::simd::simd_traits<builtin_simd_t>::scalar_t
 * \param  length     same as seqan3::simd::simd_traits<builtin_simd_t>::length
 * \param  max_length same as seqan3::simd::simd_traits<builtin_simd_t>::max_length
 * \param  simd_t     the type of the vector extension,
 * e.g. `scalar_t __attribute__ ((__vector_size__(max_length)))`
 */
#define SEQAN3_BUILTIN_SIMD_CLASS(scalar_t, length, max_length, simd_t) \
template <> \
struct builtin_simd<scalar_t, length> \
{ using type = simd_t; };

SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION(SEQAN3_BUILTIN_SIMD_CLASS)
#undef SEQAN3_BUILTIN_SIMD_CLASS

/*!\brief Helper struct for seqan3::detail::is_builtin_simd
 * \ingroup simd
 * \sa seqan3::detail::is_builtin_simd
 */
template <typename builtin_simd_t>
struct _is_builtin_simd : std::false_type
{};

/*!\brief Helper struct for seqan3::detail::is_builtin_simd
 * \ingroup simd
 * \sa seqan3::detail::is_builtin_simd
 */
template <typename builtin_simd_t>
//!\cond
    requires requires(builtin_simd_t simd){ {simd[0]} }
//!\endcond
struct _is_builtin_simd<builtin_simd_t>
{
protected:
    //!\brief seqan3::simd::simd_traits<builtin_simd_t> should be able to access *scalar_type* and *length*.
    friend class seqan3::simd::simd_traits<builtin_simd_t>;

    //!\brief The scalar type of builtin_simd_t
    using scalar_type = std::remove_reference_t<decltype(std::declval<builtin_simd_t>()[0])>;
    //!\brief The length of builtin_simd_t
    static constexpr auto length = min_viable_uint_v<sizeof(builtin_simd_t)/sizeof(scalar_type)>;

public:
    //!\brief Whether builtin_simd_t is a builtin type or not?
    //!\hideinitializer
    static constexpr bool value = std::is_same_v<builtin_simd_t, transformation_trait_or_t<builtin_simd<scalar_type, length>, void>>;
};

/*!\brief seqan3::detail::is_builtin_simd is std::true_type, **iff**
 * seqan3::detail::builtin_simd<scalar_t, length>::type is a builtin simd type.
 * \ingroup simd
 * \tparam scalar_type The underlying type of a simd vector
 * \tparam length_v The number of packed values in a simd vector
 *
 * \include test/snippet/core/simd/detail/is_builtin_simd.cpp
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html
 */
template <typename builtin_simd_t>
struct is_builtin_simd : std::bool_constant<_is_builtin_simd<builtin_simd_t>::value>
{};

/*!\brief This function specializes seqan3::detail::default_simd_max_length for
 * seqan3::detail::builtin_simd.
 * \ingroup simd
 *
 * The redefinition of *default_simd_max_length* influences the default
 * *length* (i.e., seqan3::detail::default_simd_length) of seqan3::simd::simd_type for
 * seqan3::detail::builtin_simd types.
 *
 * \attention
 * builtin_simd does not support *avx512* (64 * 8bit = 512bit) and
 * will fall back to *avx2* (32 * 8bit = 256bit).
 */
template <>
constexpr auto default_simd_max_length<builtin_simd> = []()
{
#if defined(__AVX512F__)
    return min_viable_uint_v<32u>;
#elif defined(__AVX2__)
    return min_viable_uint_v<32u>;
#elif defined(__SSE4_1__) && defined(__SSE4_2__)
    return min_viable_uint_v<16u>;
#else
    // static_assert(false, "Your code uses simd, but you didn't provide flags "
    // "(like `-march=native`, `-msse4`, `-mavx2`, etc) to your compiler "
    // "to build special binaries that make use of the special vector/simd "
    // "registers.");
    return min_viable_uint_v<0u>;
#endif
}();

} // namespace seqan3::detail

namespace seqan3::simd
{

/*!\brief This class specializes seqan3::simd::simd_traits for seqan3::detail::builtin_simd types
 * \tparam builtin_simd_t A simd type that satisfies seqan3::detail::is_builtin_simd_v<builtin_simd_t>.
 * \ingroup simd
 * \sa seqan3::simd::simd_traits for more information
 */
template <typename builtin_simd_t>
// \cond
    requires detail::is_builtin_simd<builtin_simd_t>::value
// \endcond
struct simd_traits<builtin_simd_t>
{
    //!\copydoc seqan3::simd::simd_traits::scalar_type
    using scalar_type = typename detail::_is_builtin_simd<builtin_simd_t>::scalar_type;
    //!\copydoc seqan3::simd::simd_traits::length
    static constexpr auto length = detail::_is_builtin_simd<builtin_simd_t>::length;
    //!\copydoc seqan3::simd::simd_traits::max_length
    static constexpr auto max_length = length == 1u ? length : sizeof(scalar_type) * length;

    static_assert(std::is_integral_v<scalar_type>, "For now we assume that builtin simd can only be integers");
    //!\copydoc seqan3::simd::simd_traits::mask_type
    using mask_type = typename detail::builtin_simd<std::make_signed_t<scalar_type>, length>::type;
    //!\copydoc seqan3::simd::simd_traits::swizzle_type
    using swizzle_type = typename detail::builtin_simd<uint8_t, max_length>::type;
};

} // namespace seqan3::simd
