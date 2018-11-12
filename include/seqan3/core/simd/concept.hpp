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
 * \brief Contains seqan3::simd::Simd
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

inline namespace simd
{

/*!\interface seqan3::simd::Simd <>
 * \brief The generic simd concept.
 * \ingroup simd
 *
 * \details
 *
 * seqan3::simd::Simd checks whether a given type is a simd type. One of the prerequisites is
 * that seqan3::simd::simd_traits is defined for this type.
 */
//!\cond
template <typename simd_t>
concept Simd = requires (simd_t a, simd_t b)
{
    typename simd_traits<simd_t>::scalar_type;
    typename simd_traits<simd_t>::mask_type;
    typename simd_traits<simd_t>::swizzle_type;

    // require that static member variables are defined
    requires std::Integral<decltype(simd_traits<simd_t>::length)>;
    requires std::Integral<decltype(simd_traits<simd_t>::max_length)>;

    // assume array access that returns a scalar_type type
    { a[0] } -> typename simd_traits<simd_t>::scalar_type;

    // require comparison operators
    requires std::Same<decltype(a == b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a != b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a <  b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a >  b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a <= b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a >= b), typename simd_traits<simd_t>::mask_type>;

    // require arithmetic operators
    requires std::Same<decltype(a + b), simd_t>;
    requires std::Same<decltype(a - b), simd_t>;
    requires std::Same<decltype(a * b), simd_t>;
    requires std::Same<decltype(a / b), simd_t>;
    requires std::Same<decltype(a += b), simd_t &>;
    requires std::Same<decltype(a -= b), simd_t &>;
    requires std::Same<decltype(a *= b), simd_t &>;
    requires std::Same<decltype(a /= b), simd_t &>;
};
//!\endcond

} // inline namespace simd

} // namespace seqan3
