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

#include "../alphabet.hpp"
#include "../alphabet_container.hpp"
#include "dna4.hpp"

// ------------------------------------------------------------------
// containers
// -----------------------------------------------------------------

namespace seqan3
{

using dna4_vector = std::vector<dna4>;

using dna4_string = std::basic_string<dna4, std::char_traits<dna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// -----------------------------------------------------------------

namespace seqan3::literal
{

inline dna4_vector operator "" _dna4(const char * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return dna4{dna4::char_to_value[c]};
    });

    return r;
}

inline dna4_string operator "" _dna4s(const char * s, std::size_t n)
{
    dna4_string r;
    r.resize(n);

    std::transform(s, s + n, r.begin(), [] (const char & c)
    {
        return dna4{dna4::char_to_value[c]};
    });

    return r;
}

} // namespace seqan3::literal

