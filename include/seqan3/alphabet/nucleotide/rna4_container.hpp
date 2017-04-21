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
// Author: Chenxu Pan <chenxu.pan@fu-berlin.de>
// ============================================================================

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include "../alphabet.hpp"
#include "../alphabet_container.hpp"
#include "rna4.hpp"


/*! Containers of @link rna4 @endlink
 */

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using rna4_vector = std::vector<rna4>;


/*! std::basic_string of rna4
 *
 * **NOTE** that we recommend using @link rna4_vector @endlink in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using rna4_string = std::basic_string<rna4, std::char_traits<rna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

/*! rna4 literal (returns @link rna4_vector @endlink)
 *
 * You can use this string literal to easily assign to rna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna4_vector foo{"ACGUUA"};
 *     // rna4_vector bar = "ACGUUA";
 *
 *     // but these do:
 *     rna4_vector foo{"ACGUUA"_rna4};
 *     rna4_vector bar = "ACGUUA"_rna4;
 *     auto bax = "ACGUUA"_rna4;
 *~~~~~~~~~~~~~~~
 */

inline rna4_vector operator "" _rna4(const char * s, std::size_t n)
{
    rna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].from_char(s[i]);

    return r;
}

/*! rna4 string literal (returns @link rna4_string @endlink)
 *
 * You can use this string literal to easily assign to rna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna4_string foo{"ACGUUA"};
 *     // rna4_string bar = "ACGUUA";
 *
 *     // but these do:
 *     rna4_string foo{"ACGUUA"_rna4s};
 *     rna4_string bar = "ACGUUA"_rna4s;
 *     auto bax = "ACGUUA"_rna4s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of @link rna4_string @endlink and consider using the `_rna4' literal instead.
 */

inline rna4_string operator "" _rna4s(const char * s, std::size_t n)
{
    rna4_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].from_char(s[i]);

    return r;
}

} // namespace seqan3::literal
