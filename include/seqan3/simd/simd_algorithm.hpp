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
 * \brief Contains algorithms to modify seqan3::simd.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/simd/concept.hpp>
#include <seqan3/simd/simd_traits.hpp>

#include <utility>

namespace seqan3::detail
{

//!\brief Helper function for seqan3::simd_fill.
//!\ingroup simd
template <simd_concept simd_t, size_t... I>
constexpr simd_t simd_fill_impl(typename simd_traits<simd_t>::scalar_type const scalar, std::index_sequence<I...>)
{
    return simd_t{((void)I, scalar)...};
}

//!\brief Helper function for seqan3::simd_iota.
//!\ingroup simd
template <simd_concept simd_t, typename scalar_t, scalar_t... I>
constexpr simd_t simd_iota_impl(scalar_t const offset, std::integer_sequence<scalar_t, I...>)
{
    return simd_t{(offset + I)...};
}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Fills a seqan3::simd vector with a scalar value.
 * \tparam    simd_t The simd type which satisfies seqan3::simd_concept.
 * \param[in] scalar The scalar value to fill the seqan3::simd vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/simd/simd_fill.cpp
 */
template <simd_concept simd_t>
constexpr simd_t simd_fill(typename simd_traits<simd_t>::scalar_type const scalar)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    return detail::simd_fill_impl<simd_t>(scalar, std::make_index_sequence<length>{});
}

/*!\brief Fills a seqan3::simd vector with the scalar values offset, offset+1, offset+2, ...
 * \tparam    simd_t The simd type which satisfies seqan3::simd_concept.
 * \param[in] offset The scalar offset to fill the seqan3::simd vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/simd/simd_iota.cpp
 */
template <simd_concept simd_t>
constexpr simd_t simd_iota(typename simd_traits<simd_t>::scalar_type const offset)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    using scalar_type = typename simd_traits<simd_t>::scalar_type;
    return detail::simd_iota_impl<simd_t>(offset, std::make_integer_sequence<scalar_type, length>{});
}

} // namespace seqan3
