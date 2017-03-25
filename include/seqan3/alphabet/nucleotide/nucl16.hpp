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

#include <cassert>

#include "../alphabet.hpp"

/*! The sixteen letter nucleotide alphabet
 * \ingroup alphabet
 */

namespace seqan3
{
/*! The sixteen letter nucleotide alphabet of A, B, C, D, G, H, K, M, N, R, S, T, U, V, W, Y
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 * ~~~~~~~~~~~~~~~{.cpp}
 *     nucl16 my_letter{nucl16::A};
 *     // doesn't work:
 *     // nucl16 my_letter{'A'};
 *
 *     my_letter.from_char('C'); // <- this does!
 *
 *     my_letter.from_char('E'); // converted to N internally
 *     if (my_letter.to_char() == 'N')
 *        std::cout << "yeah\n"; // "yeah";
 * ~~~~~~~~~~~~~~~
 */

struct nucl16
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
        B,
        C,
        D,
        G,
        H,
        K,
        M,
        N,
        R,
        S,
        T,
        U,
        V,
        W,
        Y,
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
    static const nucl16 A;
    static const nucl16 B;
    static const nucl16 C;
    static const nucl16 D;
    static const nucl16 G;
    static const nucl16 H;
    static const nucl16 K;
    static const nucl16 M;
    static const nucl16 N;
    static const nucl16 R;
    static const nucl16 S;
    static const nucl16 T;
    static const nucl16 U;
    static const nucl16 V;
    static const nucl16 W;
    static const nucl16 Y;
    static const nucl16 UNKNOWN;
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
    constexpr nucl16 from_char(char_type const c)
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
    constexpr nucl16 from_integral(integral_type const c)
    {
        assert(c < value_size);
        value = static_cast<internal_type>(c);
        return *this;
    }

    //! The size of the alphabet, i.e. the number of different values it can take.
    static constexpr integral_type value_size{16};

    //! @name comparison operators
    //!@{
    constexpr bool operator==(nucl16 const & rhs) const
    {
        return value == rhs.value;
    }

    constexpr bool operator!=(nucl16 const & rhs) const
    {
        return value != rhs.value;
    }

    constexpr bool operator<(nucl16 const & rhs) const
    {
        return value < rhs.value;
    }

    constexpr bool operator>(nucl16 const & rhs) const
    {
        return value > rhs.value;
    }

    constexpr bool operator<=(nucl16 const & rhs) const
    {
        return value <= rhs.value;
    }

    constexpr bool operator>=(nucl16 const & rhs) const
    {
        return value >= rhs.value;
    }
    //!@}

    //! \privatesection
    // conversion tables
    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'B',
        'C',
        'D',
        'G',
        'H',
        'K',
        'M',
        'N',
        'R',
        'S',
        'T',
        'U',
        'V',
        'W',
        'Y'
    };

    static constexpr internal_type char_to_value[256]
    {
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //                      A,                      B,                      C,
        internal_type::UNKNOWN, internal_type::A,       internal_type::B,       internal_type::C,
        //D,                    E,                      F,                      G,
        internal_type::D,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::G,
        //H,                    I,                      J,                      K,
        internal_type::H,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::K,
        //L,                    M,                      N,                      O,
        internal_type::UNKNOWN, internal_type::M,       internal_type::N,       internal_type::UNKNOWN,
        //P,                    Q,                      R,                      S,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::R,       internal_type::S,
        //T,                    U,                      V,                      W,
        internal_type::T,       internal_type::U,       internal_type::V,       internal_type::W,
        //X,                    Y,                      Z
        internal_type::UNKNOWN, internal_type::Y,       internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        //                      a,                      b,                      c,
        internal_type::UNKNOWN, internal_type::A,       internal_type::B,       internal_type::C,
        //d,                    e,                      f,                      g,
        internal_type::D,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::G,
        //h,                    i,                      j,                      k,
        internal_type::H,       internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::K,
        //l,                    m,                      n,                      o,
        internal_type::UNKNOWN, internal_type::M,       internal_type::N,       internal_type::UNKNOWN,
        //p,                    q,                      r,                      s,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::R,       internal_type::S,
        //t,                    u,                      v,                      w,
        internal_type::T,       internal_type::U,       internal_type::V,       internal_type::W,
        //x,                    y,                      z
        internal_type::UNKNOWN, internal_type::Y,       internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN
    };

};

constexpr nucl16 nucl16::A{internal_type::A};
constexpr nucl16 nucl16::B{internal_type::B};
constexpr nucl16 nucl16::C{internal_type::C};
constexpr nucl16 nucl16::D{internal_type::D};
constexpr nucl16 nucl16::G{internal_type::G};
constexpr nucl16 nucl16::H{internal_type::H};
constexpr nucl16 nucl16::K{internal_type::K};
constexpr nucl16 nucl16::M{internal_type::M};
constexpr nucl16 nucl16::N{internal_type::N};
constexpr nucl16 nucl16::R{internal_type::R};
constexpr nucl16 nucl16::S{internal_type::S};
constexpr nucl16 nucl16::T{internal_type::T};
constexpr nucl16 nucl16::U{internal_type::U};
constexpr nucl16 nucl16::V{internal_type::V};
constexpr nucl16 nucl16::W{internal_type::W};
constexpr nucl16 nucl16::Y{internal_type::Y};
constexpr nucl16 nucl16::UNKNOWN{nucl16::N};

#ifndef NDEBUG
static_assert(alphabet_concept<nucl16>);
static_assert(detail::internal_alphabet_concept<nucl16>);
#endif

}
