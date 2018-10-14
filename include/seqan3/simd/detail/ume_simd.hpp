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
 * \brief Provides seqan3::detail::ume_simd, seqan3::detail::is_ume_simd and seqan3::simd_traits<ume_simd_t>
 */

#pragma once

#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/simd/detail/default_simd_max_length.hpp>
#include <seqan3/simd/simd_traits.hpp>

#if __has_include(<umesimd/UMESimd.h>)
#include <umesimd/UMESimd.h>
#endif

//!\cond
// forward declare UME::SIMD structs
namespace UME::SIMD
{
template<typename simd_t>
class SIMDTraits;

template<typename scalar_t, uint32_t length>
class BaseVectorType;
} // namespace UME::SIMD
//!\endcond

namespace seqan3::detail
{

/*!\brief seqan3::detail::ume_simd is the class that holds the type of a
 * simd implementation called [UME::SIMD::SIMDVec]
 * (https://github.com/edanor/umesimd).
 * \ingroup simd
 * \tparam scalar_type The underlying type of a simd vector
 * \tparam length_v The number of packed values in a simd vector
 *
 * \include test/snippet/simd/detail/ume_simd.cpp
 * \attention This class itself only delegates to a [UME::SIMD::SIMDVec]
 * (https://github.com/edanor/umesimd) type.
 * \sa https://github.com/edanor/umesimd
 */
template <typename scalar_t, size_t length>
struct ume_simd;

/*!\ingroup simd
 * \copydoc seqan3::detail::ume_simd
 * \sa seqan3::detail::ume_simd
 */
template <typename scalar_t, size_t length>
//!\cond
    requires requires() { typename UME::SIMD::BaseVectorType<scalar_t, length>::BASE_T; }
//!\endcond
struct ume_simd<scalar_t, length>
{
    //!\brief The delegated *UME::SIMD* type.
    using type = typename UME::SIMD::BaseVectorType<scalar_t, length>::BASE_T;
};

/*!\brief Helper concept for seqan3::detail::is_ume_simd
 * \ingroup simd
 * \sa seqan3::detail::is_ume_simd
 */
template <typename ume_simd_t>
concept _is_ume_simd = requires()
{
    typename UME::SIMD::SIMDTraits<ume_simd_t>::SCALAR_T;
};

/*!\brief seqan3::detail::is_ume_simd is std::true_type, **iff**
 * seqan3::detail::ume_simd<scalar_t,length>::type is a [UME::SIMD::SIMDVec]
 * (https://github.com/edanor/umesimd) type.
 * \ingroup simd
 * \tparam scalar_type The underlying type of a simd vector
 * \tparam length_v The number of packed values in a simd vector
 *
 * \include test/snippet/simd/detail/is_ume_simd.cpp
 * \sa https://github.com/edanor/umesimd
 */
template <typename ume_simd_t, typename = std::void_t<> >
struct is_ume_simd : std::bool_constant<_is_ume_simd<ume_simd_t>>
{};

/*!\brief This function specializes seqan3::detail::default_simd_max_length for
 * seqan3::detail::ume_simd.
 * \ingroup simd
 *
 * The redefinition of seqan3::detail::default_simd_max_length influences the default
 * *length* (i.e., seqan3::detail::default_simd_length) of seqan3::simd for
 * seqan3::detail::ume_simd types.
 */
template <>
constexpr auto default_simd_max_length<ume_simd> = []()
{
#if defined(__AVX512F__)
    return min_viable_uint_v<64u>;
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

namespace seqan3
{

/*!\brief This class specializes seqan3::simd_traits for seqan3::detail::ume_simd types
 * \tparam ume_simd_t A simd type that satisfies seqan3::detail::is_ume_simd<ume_simd_t>.
 * \ingroup simd
 * \sa seqan3::simd_traits for more information
 */
template <typename ume_simd_t>
// \cond
    requires detail::is_ume_simd<ume_simd_t>::value
// \endcond
struct simd_traits<ume_simd_t>
{
    //!\copydoc seqan3::simd_traits::scalar_type
    using scalar_type = typename UME::SIMD::SIMDTraits<ume_simd_t>::SCALAR_T;
    //!\copydoc seqan3::simd_traits::length
    static constexpr auto length = detail::min_viable_uint_v<ume_simd_t::length()>;
    //!\copydoc seqan3::simd_traits::max_length
    static constexpr auto max_length = sizeof(scalar_type) * length;
    //!\copydoc seqan3::simd_traits::mask_type
    using mask_type = typename UME::SIMD::SIMDTraits<ume_simd_t>::MASK_T;
    //!\copydoc seqan3::simd_traits::swizzle_type
    using swizzle_type = typename UME::SIMD::SIMDTraits<ume_simd_t>::SWIZZLE_T;
};

} // namespace seqan3
