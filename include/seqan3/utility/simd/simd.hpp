// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::simd::simd_type
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/simd/detail/default_simd_backend.hpp>

namespace seqan3
{

inline namespace simd
{

/*!\brief seqan3::simd::simd_type encapsulates simd vector types, which can be manipulated
 * by simd operations.
 * \ingroup utility_simd
 * \tparam scalar_t The underlying type of a simd vector
 * \tparam length The number of packed values in a simd vector
 * \tparam simd_backend The simd backend to use, e.g.
 * seqan3::detail::builtin_simd
 *
 * \include test/snippet/utility/simd/simd.cpp
 * \attention
 * seqan3::simd::simd_type may not support *float* types depending on the selected backend.
 *
 * All implementations support *[u]intX_t* types, e.g. *uint8_t*.
 *
 * ###Helper types
 *   seqan3::simd::simd_type_t as a shorthand for seqan3::simd::simd_type::type
 * \sa https://en.wikipedia.org/wiki/SIMD What is SIMD conceptually?
 * \sa https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions Which SIMD architectures exist?
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html Underlying technique of *seqan3::detail::builtin_simd types*.
 * \sa https://github.com/edanor/umesimd Underlying library of *seqan3::detail::ume_simd* types.
 * \sa https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html Instruction sets and their low-level intrinsics.
 */
template <typename scalar_t,
          size_t length = detail::default_simd_length<scalar_t, detail::default_simd_backend>,
          template <typename scalar_t_, size_t length_> typename simd_backend = detail::default_simd_backend>
struct simd_type : simd_backend<scalar_t, length>
{
    //!\brief The actual simd type.
    using type = typename simd_backend<scalar_t, length>::type;
};

//!\brief Helper type of seqan3::simd::simd_type
//!\ingroup utility_simd
template <typename scalar_t,
          size_t length = detail::default_simd_length<scalar_t, detail::default_simd_backend>,
          template <typename scalar_t_, size_t length_> typename simd_backend = detail::default_simd_backend>
using simd_type_t = typename simd_type<scalar_t, length, simd_backend>::type;

} // namespace simd

} // namespace seqan3
