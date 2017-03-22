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
// Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================

#pragma once

#include "../alphabet.hpp"

#include <tuple>

// assume sequential mapping w.r.t. ascii alphabet, only (char_start, val_start) and (char_end, val_end) differ
namespace seqan3
{
    // integer score in [0 .. 41], char score in  ['!' .. 'J']
    struct illumina18
    {
        using char_type = char;
        using integral_type = uint8_t;
        using phred_type = int8_t;

        // todo: make this offsets concept
        static constexpr char_type offset_char{'!'};
        static constexpr integral_type offset_phred{0};

        // the value
        integral_type value;

        // implicit compatibility to inner_type
        constexpr illumina18 & operator =(integral_type const c)
        {
            value = c;
        }

        // explicit compatibility to char
        explicit constexpr operator char() const
        {
            return to_char();
        }

        constexpr char_type to_char() const
        {
            return value + offset_char;
        }

        constexpr illumina18 from_char(char_type const c)
        {
            value = c - '!';
            return *this;
        }

        // explicit compatibility to integral
        constexpr integral_type to_integral() const
        {
            return value;
        }

        constexpr illumina18 from_integral(uint8_t const c)
        {
            value = c;
            return *this;
        }

        constexpr illumina18 from_phred(int8_t const p)
        {
            value = p - offset_phred;
            return *this;
        }

        constexpr phred_type to_phred(uint8_t const p)
        {
            return value + offset_phred;
        }

        //
        static constexpr uint8_t value_size{42};

    };

}
