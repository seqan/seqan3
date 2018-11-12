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
 * \brief Provides seqan3::debug_stream overload for seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

namespace seqan3
{

/*!\brief Overload for debug_stream for simd types.
 * \ingroup simd
 */
template <typename simd_t>
//!\cond
    requires simd::Simd<remove_cvref_t<simd_t>>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & s, simd_t && simd)
{
    using simd_type = remove_cvref_t<simd_t>;
    constexpr size_t length = simd::simd_traits<simd_type>::length;
    using scalar_type = typename simd::simd_traits<simd_type>::scalar_type;

    std::array<scalar_type, length> array{};
    for (size_t i = 0; i < length; ++i)
        array[i] = simd[i];
    s << array;
    return s;
}

} // namespace seqan3
