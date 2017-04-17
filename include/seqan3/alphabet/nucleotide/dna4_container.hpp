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
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/range.hpp>
#include "dna4.hpp"


/*! Containers of @link dna4 @endlink
 */

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using dna4_vector = std::vector<dna4>;


/*! std::basic_string of dna4
 *
 * **NOTE** that we recommend using @link dna4_vector @endlink in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using dna4_string = std::basic_string<dna4, std::char_traits<dna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

/*! dna4 literal (returns @link dna4_vector @endlink)
 *
 * You can use this string literal to easily assign to dna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna4_vector foo{"ACGTTA"};
 *     // dna4_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     dna4_vector foo{"ACGTTA"_dna4};
 *     dna4_vector bar = "ACGTTA"_dna4;
 *     auto bax = "ACGTTA"_dna4;
 *~~~~~~~~~~~~~~~
 */

inline dna4_vector operator "" _dna4(const char * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].from_char(s[i]);

    return r;
}

/*! dna4 string literal (returns @link dna4_string @endlink)
 *
 * You can use this string literal to easily assign to dna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna4_string foo{"ACGTTA"};
 *     // dna4_string bar = "ACGTTA";
 *
 *     // but these do:
 *     dna4_string foo{"ACGTTA"_dna4s};
 *     dna4_string bar = "ACGTTA"_dna4s;
 *     auto bax = "ACGTTA"_dna4s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of @link dna4_string @endlink and consider using the `_dna4' literal instead.
 */

inline dna4_string operator "" _dna4s(const char * s, std::size_t n)
{
    dna4_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].from_char(s[i]);

    return r;
}

} // namespace seqan3::literal
