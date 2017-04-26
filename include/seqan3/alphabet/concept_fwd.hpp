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

/*!\file alphabet/concept_fwd.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Core alphabet concept metafunction base classes.
 *
 * Include this file, if you implement an alphabet type with Free/global function
 * and metafunction interface.
 */

#pragma once

#include <iostream>
#include <string>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

// ============================================================================
// Free/Global interface wrappers
// ============================================================================

//!\addtogroup alphabet
//!\{

// ------------------------------------------------------------------
// type metafunctions
// ------------------------------------------------------------------

//!\brief Type metafunction that shall return the type of an alphabet in char representation. [base type]
//!\tparam alphabet_type The unspecialized and unconstrained template takes any type.
template <typename alphabet_type>
struct underlying_char;

//!\brief Shortcut for seqan3::underlying_char
template <typename alphabet_type>
//!\relates seqan3::underlying_char
using underlying_char_t = typename underlying_char<alphabet_type>::type;

//!\brief Type metafunction that shall return the type of an alphabet in rank representation. [base type]
//!\tparam alphabet_type Must provide a `rank_type` member type.
template <typename alphabet_type>
struct underlying_rank;

//!\brief Shortcut for seqan3::underlying_rank
//!\relates seqan3::underlying_rank
template <typename alphabet_type>
using underlying_rank_t = typename underlying_rank<alphabet_type>::type;

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

//!\brief Value metafunction that shall return the size of an alphabet. [base type]
//!\tparam alphabet_type The unspecialized and unconstrained template takes any type.
template <typename alphabet_type>
struct alphabet_size;

//!\brief Shortcut for seqan3::alphabet_size
//!\relates seqan3::alphabet_size
template <typename alphabet_type>
constexpr auto alphabet_size_v = alphabet_size<alphabet_type>::value;

//!\}

} // namespace seqan3
