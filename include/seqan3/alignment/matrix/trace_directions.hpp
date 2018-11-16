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
 * \brief Contains the declaration of seqan3::detail::trace_directions.
 */

#pragma once

#include <seqan3/core/add_enum_bitwise_operators.hpp>

namespace seqan3::detail
{

/*!\brief The possible directions a trace can have. The values can be combined by the logical `|`-operator.
 * \ingroup alignment_matrix
 * \sa seqan3::add_enum_bitwise_operators <seqan3::detail::trace_directions> enables combining enum values.
 * \sa seqan3::detail::alignment_trace_matrix implementations use this enum as matrix entry type.
 */
enum struct trace_directions : uint8_t
{
    //!\brief No trace
    none      = 0b0000,
    //!\brief Trace comes from the diagonal entry.
    diagonal  = 0b0001,
    //!\brief Trace comes from the above entry.
    up        = 0b0010,
    //!\brief Trace comes from the left entry.
    left      = 0b0100
};

} // namespace seqan3::detail

namespace seqan3
{
//!\brief Enable bitwise operators for enum seqan3::detail::trace_directions.
//!\ingroup alignment_matrix
template <>
constexpr bool add_enum_bitwise_operators<seqan3::detail::trace_directions> = true;
} // namespace seqan3
