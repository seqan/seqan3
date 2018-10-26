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
 * \brief Contains seqan3::simd::simd_type
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/simd/detail/default_simd_backend.hpp>
#include <seqan3/core/simd/detail/default_simd_length.hpp>

namespace seqan3
{

inline namespace simd
{

/*!\brief seqan3::simd::simd_type encapsulates simd vector types, which can be manipulated
 * by simd operations.
 * \ingroup simd
 * \tparam scalar_t The underlying type of a simd vector
 * \tparam length The number of packed values in a simd vector
 * \tparam simd_backend The simd backend to use, e.g.
 * seqan3::detail::builtin_simd
 *
 * \include test/snippet/core/simd/simd.cpp
 * \attention
 * seqan3::simd::simd_type may not support *float* types depending on the selected backend.
 *
 * All implementations support *[u]intX_t* types, e.g. *uint8_t*.
 *
 * \par Helper types
 *   seqan3::simd::simd_type_t as a shorthand for seqan3::simd::simd_type::type
 * \sa https://en.wikipedia.org/wiki/SIMD What is SIMD conceptually?
 * \sa https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions Which SIMD architectures exist?
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html Underlying technique of *seqan3::detail::builtin_simd types*.
 * \sa https://github.com/edanor/umesimd Underlying library of *seqan3::detail::ume_simd* types.
 * \sa https://software.intel.com/sites/landingpage/IntrinsicsGuide Instruction sets and their low-level intrinsics.
 */
template <typename scalar_t,
          size_t length = detail::default_simd_length<scalar_t, detail::default_simd_backend>,
          typename simd_backend = detail::default_simd_backend<scalar_t, length>>
struct simd_type : simd_backend
{
    //!\brief The actual simd type.
    using type = typename simd_backend::type;
};

//!\brief Helper type of seqan3::simd::simd_type
//!\ingroup simd
template <typename scalar_t,
          size_t length = detail::default_simd_length<scalar_t, detail::default_simd_backend>,
          typename simd_backend = detail::default_simd_backend<scalar_t, length>>
using simd_type_t = typename simd_type<scalar_t, length, simd_backend>::type;

} // inline namespace simd

} // namespace seqan3
