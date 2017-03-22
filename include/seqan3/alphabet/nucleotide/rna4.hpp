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

#include "../alphabet.hpp"

namespace seqan3
{

struct rna4
{
    using char_type = char;
    using integral_type = uint8_t;

    // strictly typed enum, unfortunately with scope
    enum struct c_type : integral_type
    {
        A,
        C,
        G,
        U,
        UNKNOWN = A
    };

    // implicit compatibility to inner_type
    constexpr rna4 & operator =(c_type const c)
    {
        value = c;
    }

    constexpr rna4 from_char(char const c)
    {
        value = char_to_value[c];
        return *this;
    }

    constexpr rna4 from_integral(integral_type const c)
    {
        value = static_cast<c_type>(c);
        return *this;
    }

    static constexpr char value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'U'
    };
    
};

// shall fulfill Alphabet concept
static_assert(alphabet_concept<rna4>);
static_assert(rna4{dna4::A} == dna4{});
static_assert(rna4{dna4::A} == dna4::A);
// static_assert(rna4{'A'} == 'A');
static_assert(static_cast<char>(rna4{dna4::C}) == 'C');
static_assert(rna4{dna4::A} < dna4{dna4::C});

}
