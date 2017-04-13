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
#include "nucl16.hpp"

/*! Containers of @link nucl16 @endlink
 */

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using nucl16_vector = std::vector<nucl16>;

/*! std::basic_string of nucl16
 *
 * **NOTE** that we recommend using @link nucl16_vector @endlink in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using nucl16_string = std::basic_string<nucl16, std::char_traits<nucl16>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

/*! nucl16 literal (returns @link nucl16_vector @endlink)
 *
 * You can use this string literal to easily assign to nucl16_vector:
 *
 *     // these don't work:
 *     // nucl16_vector foo{"ACGTTA"};
 *     // nucl16_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     nucl16_vector foo{"ACGTTA"_nucl16};
 *     nucl16_vector bar = "ACGTTA"_nucl16;
 *     auto bax = "ACGTTA"_nucl16;
 */

inline nucl16_vector operator "" _nucl16(const char * s, std::size_t n)
{
    nucl16_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].from_char(s[i]);

    return r;
}

/*! nucl16 string literal (returns @link nucl16_string @endlink)
 *
 * You can use this string literal to easily assign to nucl16_vector:
 *
 *     // these don't work:
 *     // nucl16_string foo{"ACGTTA"};
 *     // nucl16_string bar = "ACGTTA";
 *
 *     // but these do:
 *     nucl16_string foo{"ACGTTA"_nucl16s};
 *     nucl16_string bar = "ACGTTA"_nucl16s;
 *     auto bax = "ACGTTA"_nucl16s;
 *
 * Please note the limitations of @link nucl16_string @endlink and consider using the `_nucl16' literal instead.
 */

inline nucl16_string operator "" _nucl16s(const char * s, std::size_t n)
{
    nucl16_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].from_char(s[i]);

    return r;
}
} // namespace seqan3::literal

