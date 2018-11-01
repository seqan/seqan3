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
 * \brief Contains seqan3::simd::simd_traits
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

inline namespace simd
{

/*!\brief seqan3::simd::simd_traits is the trait class that provides uniform interface
 * to the properties of simd_t types.
 * \ingroup simd
 * \tparam simd_t The simd type that satisfies seqan3::simd::simd_concept.
 *
 * The class defines the following member variables and types:
 * * scalar_type - the underlying type of a simd vector
 * * length - the number of packed values in a simd vector
 * * max_length - the maximum number of packable values in a simd vector, if the underlying type would be [u]int8_t
 * * mask_type - the type returned by comparison operators
 * * swizzle_type - the type used to define how to swizzle a simd vector
 *
 * \include test/snippet/core/simd/simd_traits.cpp
 */
template <typename simd_t>
struct simd_traits
#if SEQAN3_DOXYGEN_ONLY(1)0
{
    /*!\brief The underlying type of a simd vector (is not defined if *simd_t*
     * does not satisfy *seqan3::simd::simd_concept*)
     */
    using scalar_type = IMPLEMENTATION_DEFINED;
    /*!\brief The number of packed values in a simd vector (is not defined if
     * *simd_t* does not satisfy *seqan3::simd::simd_concept*)
     */
    static constexpr auto length = IMPLEMENTATION_DEFINED;
    /*!\brief The maximum number of packable values in a simd vector, if the
     * underlying type would be *[u]int8_t* (is not defined if *simd_t* does not
     * satisfy *seqan3::simd::simd_concept*)
     */
    static constexpr auto max_length = IMPLEMENTATION_DEFINED;
    /*!\brief The type returned by comparison operators (is not defined if
     * *simd_t* does not satisfy *seqan3::simd::simd_concept*)
     */
    using mask_type = IMPLEMENTATION_DEFINED;
    /*!\brief The type used to define how to swizzle a simd vector (is not
     * defined if *simd_t* does not satisfy *seqan3::simd::simd_concept*)
     */
    using swizzle_type = IMPLEMENTATION_DEFINED;
}
#endif
;

} // inline namespace simd

} // namespace seqan3
