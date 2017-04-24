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
// Author: Sara Hetzel <sara.hetzel AT fu-berlin.de>
// ============================================================================

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/range.hpp>
#include "aa27.hpp"

/*! Containers of @link aa27 @endlink
 */

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using aa27_vector = std::vector<aa27>;

/*! std::basic_string of aa27
 *
 * **NOTE** that we recommend using @link aa27_vector @endlink in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using aa27_string = std::basic_string<aa27, std::char_traits<aa27>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

/*! aa27 literal (returns @link aa27_vector @endlink)
 *
 * You can use this string literal to easily assign to aa27_vector:
 *
 *     // these don't work:
 *     // aa27_vector foo{"ABFUYR"};
 *     // aa27_vector bar = "ABFUYR";
 *
 *     // but these do:
 *     aa27_vector foo{"ABFUYR"_aa27};
 *     aa27_vector bar = "ABFUYR"_aa27;
 *     auto bax = "ABFUYR"_aa27;
 */

inline aa27_vector operator "" _aa27(const char * s, std::size_t n)
{
    aa27_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*! aa27 string literal (returns @link aa27_string @endlink)
 *
 * You can use this string literal to easily assign to aa27_vector:
 *
 *     // these don't work:
 *     // aa27_string foo{"ABFUYR"};
 *     // aa27_string bar = "ABFUYR";
 *
 *     // but these do:
 *     aa27_string foo{"ABFUYR"_aa27s};
 *     aa27_string bar = "ABFUYR"_aa27s;
 *     auto bax = "ABFUYR"_aa27s;
 *
 * Please note the limitations of @link aa27_string @endlink and consider using the `_aa27' literal instead.
 */

inline aa27_string operator "" _aa27s(const char * s, std::size_t n)
{
    aa27_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal

