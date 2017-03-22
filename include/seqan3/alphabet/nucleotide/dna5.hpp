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
// Implementation of the biological dna5 alphabet.
// ==========================================================================

#pragma once

#include "../alphabet.hpp"

namespace seqan3
{
struct dna5
{
    using char_type = char;
    using integral_type = uint8_t;

    // strictly typed enum, unfortunately with scope
    enum struct c_type : integral_type
    {
        A,
        C,
        G,
        T,
        N,
        U = T,
        UNKNOWN = N
    };

    // import into local scope
    static constexpr c_type A{c_type::A};
    static constexpr c_type C{c_type::C};
    static constexpr c_type G{c_type::G};
    static constexpr c_type T{c_type::T};
    static constexpr c_type N{c_type::N};
    static constexpr c_type U{c_type::U};
    static constexpr c_type UNKNOWN{c_type::UNKNOWN};

    // the value
    c_type value;

    // implicit compatibility to inner_type
    constexpr dna5 & operator =(c_type const c)
    {
        value = c;
    }
    constexpr operator c_type() const
    {
        return value;
    }

    // explicit compatibility to char
    explicit constexpr operator char_type() const
    {
        return to_char();
    }
    constexpr char_type to_char() const
    {
        return value_to_char[static_cast<integral_type>(value)];
    }

    constexpr dna5 from_char(char_type const c)
    {
        value = char_to_value[c];
        return *this;
    }

    // explicit compatibility to integral
    constexpr integral_type to_integral() const
    {
        return static_cast<integral_type>(value);
    }

    constexpr dna5 from_integral(integral_type const c)
    {
        value = static_cast<c_type>(c);
        return *this;
    }

    // conversion tables
    static constexpr integral_type value_size{5};

    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'T',
        'N'
    };

    static constexpr c_type char_to_value[256]
    {
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//5
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//30
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//55
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        //A,             B,               C,               D,               E,
        c_type::A,       c_type::UNKNOWN, c_type::C,       c_type::UNKNOWN, c_type::UNKNOWN,
        //F,             G,               H,               I,               J,
        c_type::UNKNOWN, c_type::G,       c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        //K,             L,               M,               N,               O,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::N,       c_type::UNKNOWN,//80
        // P,            Q,               R,               S,               T,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::T,
        //U,             V,               W,               X,               Y,
        c_type::T,       c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        //Z
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        //                                a,               b,               c,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::A,       c_type::UNKNOWN, c_type::C,
        //d,             e,               f,               g,               h,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::G,       c_type::UNKNOWN,//105
        //i,             j,               k,               l,               m,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        //n,             o,               p,               q,               r,
        c_type::N,       c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        //s,             t,               u,               v,               w,
        c_type::UNKNOWN, c_type::T,       c_type::T,       c_type::UNKNOWN, c_type::UNKNOWN,
        //x,             y,               z
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//130
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//155
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//180
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//205
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//230
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,
        c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN, c_type::UNKNOWN,//255
        c_type::UNKNOWN
    };
};
} // namespace seqan3
