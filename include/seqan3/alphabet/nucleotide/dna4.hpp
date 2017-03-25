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

#include <cassert>

#include "../alphabet.hpp"

/*! The four letter DNA alphabet
 * \ingroup alphabet
 */

namespace seqan3
{

/*! The four letter DNA alphabet of A,C,G,T
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 *     dna4 my_letter{dna4::A};
 *     // doesn't work:
 *     // dna4 my_letter{'A'};
 *
 *     my_letter.from_char('C'); // <- this does!
 *
 *     my_letter.from_char('F'); // converted to A internally
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 */

struct dna4
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
        U = T,
        UNKNOWN = A
    };
    internal_type value;
    //! \publicsection

    // import internal_types values into local scope:

    /*! @name letter values
     * Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * *Don't worry about the `internal_type`.*
     */
    //!@{
    static const dna4 A;
    static const dna4 C;
    static const dna4 G;
    static const dna4 T;
    static const dna4 U;
    static const dna4 UNKNOWN;
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
    constexpr dna4 from_char(char_type const c)
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
    constexpr dna4 from_integral(integral_type const c)
    {
        assert(c < value_size);
        value = static_cast<internal_type>(c);
        return *this;
    }

    //! The size of the alphabet, i.e. the number of different values it can take.
    static constexpr integral_type value_size{4};

    //! @name comparison operators
    //!@{
    constexpr bool operator==(dna4 const & rhs) const
    {
        return value == rhs.value;
    }

    constexpr bool operator!=(dna4 const & rhs) const
    {
        return value != rhs.value;
    }

    constexpr bool operator<(dna4 const & rhs) const
    {
        return value < rhs.value;
    }

    constexpr bool operator>(dna4 const & rhs) const
    {
        return value > rhs.value;
    }

    constexpr bool operator<=(dna4 const & rhs) const
    {
        return value <= rhs.value;
    }

    constexpr bool operator>=(dna4 const & rhs) const
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
        'T'
    };

    static constexpr internal_type char_to_value[256]
    {
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //15
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //31
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //47
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //63
        //              ,                      A,                      B,                      C,
        internal_type::UNKNOWN, internal_type::A,       internal_type::UNKNOWN, internal_type::C,
        //             D,                      E,                      F,                      G,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::G,
        //             H,                      I,                      J,                      K,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //             L,                      M,                      N,                      O,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //79
        //             P,                      Q,                      R,                      S,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //             T,                      U,                      V,                      W,
        internal_type::T,       internal_type::T,       internal_type::UNKNOWN, internal_type::UNKNOWN,
        //             X,                      Y,                      Z,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //95
        //              ,                      a,                      b,                      c,
        internal_type::UNKNOWN, internal_type::A,       internal_type::UNKNOWN, internal_type::C,
        //             d,                      e,                      f,                      g,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::G,
        //             h,                      i,                      j,                      k,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //             l,                      m,                      n,                      o,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //111
        //             p,                      q,                      r,                      s,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //             t,                      u,                      v,                      w,
        internal_type::T,       internal_type::T,       internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //127
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //143
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //159
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //175
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //191
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //207
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //223
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, //239
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN  //255
    };
};

constexpr dna4 dna4::A{internal_type::A};
constexpr dna4 dna4::C{internal_type::C};
constexpr dna4 dna4::G{internal_type::G};
constexpr dna4 dna4::T{internal_type::T};
constexpr dna4 dna4::U{dna4::T};
constexpr dna4 dna4::UNKNOWN{dna4::A};

#ifndef NDEBUG
static_assert(alphabet_concept<dna4>);
static_assert(detail::internal_alphabet_concept<dna4>);
#endif

} // namespace seqan3
