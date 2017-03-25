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

#include "../alphabet.hpp"

/*! The five letter DNA alphabet
 * \ingroup alphabet
 */

namespace seqan3
{

/*! The five letter DNA alphabet of A,C,G,T,N
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 *     dna5 my_letter{dna5::A};
 *     // doesn't work:
 *     // dna5 my_letter{'A'};
 *
 *     my_letter.from_char('C'); // <- this does!
 *
 *     my_letter.from_char('F'); // converted to N internally
 *     if (my_letter.to_char() == 'N')
 *        std::cout << "yeah\n"; // "yeah";
 */

struct dna5
{
    //! the type of the alphabet when converted to char (e.g. via @link to_char @endlink)
    using char_type = char;
    //! the type of the alphabet when represented as a number (e.g. via @link to_integral @endlink)
    using integral_type = uint8_t;

    // strictly typed enum, unfortunately with scope
    //! \privatesection
    enum struct internal_type : integral_type
    {
        A,
        C,
        G,
        T,
        N,
        U = T,
        UNKNOWN = N
    };
    internal_type value;
    //! \publicsection

    // import internal_types values into local scope:

    /*! @name letter values
     * Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * *Don't worry about the `internal_type`.*
     */
    //!@{
    static const dna5 A;
    static const dna5 C;
    static const dna5 G;
    static const dna5 T;
    static const dna5 N;
    static const dna5 U;
    static const dna5 UNKNOWN;
    //!@}

    //! ability to cast to @link char_type @endlink **explicitly**.
    explicit constexpr operator char_type() const
    {
        return to_char();
    }

    //! return the letter as a character of @link char_type @endlink.
    constexpr char_type to_char() const
    {
        return value_to_char[static_cast<integral_type>(value)];
    }

    //! assign from a character
    constexpr dna5 & from_char(char_type const c)
    {
        value = char_to_value[c];
        return *this;
    }

    //! return the letter's numeric value or rank in the alphabet
    constexpr integral_type to_integral() const
    {
        return static_cast<integral_type>(value);
    }

    //! assign from a numeric value
    constexpr dna5 & from_integral(integral_type const c)
    {
        value = static_cast<internal_type>(c);
        return *this;
    }

    //! The size of the alphabet, i.e. the number of different values it can take.
    static constexpr integral_type value_size{5};

    //! @name comparison operators
    //!@{
    constexpr bool operator==(dna5 const & rhs) const
    {
        return value == rhs.value;
    }

    constexpr bool operator!=(dna5 const & rhs) const
    {
        return value != rhs.value;
    }

    constexpr bool operator<(dna5 const & rhs) const
    {
        return value < rhs.value;
    }

    constexpr bool operator>(dna5 const & rhs) const
    {
        return value > rhs.value;
    }

    constexpr bool operator<=(dna5 const & rhs) const
    {
        return value <= rhs.value;
    }

    constexpr bool operator>=(dna5 const & rhs) const
    {
        return value >= rhs.value;
    }
    //!@}

    //! \privatesection
    // conversion tables
    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'T',
        'N'
    };

    static constexpr internal_type char_to_value[256]
    {
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//4
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//29
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//54
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //A,             B,               C,               D,               E,
        internal_type::A,       internal_type::UNKNOWN, internal_type::C,       internal_type::UNKNOWN, internal_type::UNKNOWN,
        //F,             G,               H,               I,               J,
        internal_type::UNKNOWN, internal_type::G,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //K,             L,               M,               N,               O,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::N,       internal_type::UNKNOWN,//79
        // P,            Q,               R,               S,               T,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::T,
        //U,             V,               W,               X,               Y,
        internal_type::T,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //Z
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //                                a,               b,               c,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::A,       internal_type::UNKNOWN, internal_type::C,
        //d,             e,               f,               g,               h,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::G,       internal_type::UNKNOWN,//104
        //i,             j,               k,               l,               m,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //n,             o,               p,               q,               r,
        internal_type::N,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //s,             t,               u,               v,               w,
        internal_type::UNKNOWN, internal_type::T,       internal_type::T,       internal_type::UNKNOWN, internal_type::UNKNOWN,
        //x,             y,               z
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//129
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//154
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//179
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//204
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//229
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,//254
        internal_type::UNKNOWN
    };
};

constexpr dna5 dna5::A{internal_type::A};
constexpr dna5 dna5::C{internal_type::C};
constexpr dna5 dna5::G{internal_type::G};
constexpr dna5 dna5::T{internal_type::T};
constexpr dna5 dna5::N{internal_type::N};
constexpr dna5 dna5::U{dna5::T};
constexpr dna5 dna5::UNKNOWN{dna5::N};

#ifndef NDEBUG
static_assert(alphabet_concept<dna5>);
static_assert(detail::internal_alphabet_concept<dna5>);
#endif

} // namespace seqan3
