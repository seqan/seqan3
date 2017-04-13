// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================
// Author: David Heller <david.heller@fu-berlin.de>
// ==========================================================================

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/range.hpp>
#include "dna5.hpp"

/*! Containers of @link dna5 @endlink
 */

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using dna5_vector = std::vector<dna5>;

/*! std::basic_string of dna5
 *
 * **NOTE** that we recommend using @link dna5_vector @endlink in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using dna5_string = std::basic_string<dna5, std::char_traits<dna5>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

/*! dna5 literal (returns @link dna5_vector @endlink)
 *
 * You can use this string literal to easily assign to dna5_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna5_vector foo{"ACGTTA"};
 *     // dna5_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     dna5_vector foo{"ACGTTA"_dna5};
 *     dna5_vector bar = "ACGTTA"_dna5;
 *     auto bax = "ACGTTA"_dna5;
 *~~~~~~~~~~~~~~~
 */

inline dna5_vector operator "" _dna5(const char * s, std::size_t n)
{
    dna5_vector r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return dna5{}.from_char(c);
    });

    return r;
}

/*! dna5 string literal (returns @link dna5_string @endlink)
 *
 * You can use this string literal to easily assign to dna5_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna5_string foo{"ACGTTA"};
 *     // dna5_string bar = "ACGTTA";
 *
 *     // but these do:
 *     dna5_string foo{"ACGTTA"_dna5s};
 *     dna5_string bar = "ACGTTA"_dna5s;
 *     auto bax = "ACGTTA"_dna5s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of @link dna5_string @endlink and consider using the `_dna5' literal instead.
 */

inline dna5_string operator "" _dna5s(const char * s, std::size_t n)
{
    dna5_string r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return dna5{}.from_char(c);
    });

    return r;
}

} // namespace seqan3::literal
