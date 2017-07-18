// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \ingroup structure
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::structure_concept.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// auxiliary metafunction
// ============================================================================

namespace seqan3::detail
{

//!\brief Metafunction that indicates whether an alphabet is a structure alphabet.
//!\ingroup structure
template <typename type>
struct is_structure : public std::false_type
{};

//!\brief Shortcut for seqan3::detail::is_structure.
//!\ingroup structure
template <typename type>
constexpr bool is_structure_v = is_structure<type>::value;

} // namespace seqan3::detail

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

//!\brief A concept that indicates whether an alphabet represents structure.
//!\ingroup structure
//!\details Refines the seqan3::structure_concept.
template <typename type>
concept bool structure_concept = alphabet_concept<type> && detail::is_structure_v<type>;

} // namespace seqan3
