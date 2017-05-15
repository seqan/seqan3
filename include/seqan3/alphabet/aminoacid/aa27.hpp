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

#include <seqan3/alphabet/concept.hpp>

/*! The twenty-seven letter amino acid alphabet
 * \ingroup alphabet
 */

namespace seqan3
{
/*! The twenty-seven letter amino acid alphabet A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X,
 * Y, Z, *
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 * ~~~~~~~~~~~~~~~{.cpp}
 *     aa27 my_letter{aa27::A};
 *     // doesn't work:
 *     // aa27 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('?'); // converted to X internally
 *     if (my_letter.to_char() == 'X')
 *        std::cout << "yeah\n"; // "yeah";
 * ~~~~~~~~~~~~~~~
 */

struct aa27
{
    //! the type of the alphabet when converted to char (e.g. via @link to_char @endlink)
    using char_type = char;

    //! the type of the alphabet when represented as a number (e.g. via @link to_rank @endlink)
    using rank_type = uint8_t;

    // strictly typed enum, unfortunately with scope
    //! \privatesection
    enum struct internal_type : rank_type
    {
        A,
        B,
        C,
        D,
        E,
        F,
        G,
        H,
        I,
        J,
        K,
        L,
        M,
        N,
        O,
        P,
        Q,
        R,
        S,
        T,
        U,
        V,
        W,
        X,
        Y,
        Z,
        TERMINATOR,
        UNKNOWN = X
    };
    internal_type value;
    //! \publicsection

    // import internal_types values into local scope:

    /*! @name letter values
     * Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * *Don't worry about the `internal_type`.*
     */
    //!@{

    // import into local scope
    static const aa27 A;
    static const aa27 B;
    static const aa27 C;
    static const aa27 D;
    static const aa27 E;
    static const aa27 F;
    static const aa27 G;
    static const aa27 H;
    static const aa27 I;
    static const aa27 J;
    static const aa27 K;
    static const aa27 L;
    static const aa27 M;
    static const aa27 N;
    static const aa27 O;
    static const aa27 P;
    static const aa27 Q;
    static const aa27 R;
    static const aa27 S;
    static const aa27 T;
    static const aa27 U;
    static const aa27 V;
    static const aa27 W;
    static const aa27 X;
    static const aa27 Y;
    static const aa27 Z;
    static const aa27 TERMINATOR;
    static const aa27 UNKNOWN;
    //!@}

    //! ability to cast to @link char_type @endlink **explicitly**.
    explicit constexpr operator char_type() const
    {
        return to_char();
    }

    //! return the letter as a character of @link char_type @endlink.
    constexpr char_type to_char() const
    {
        return value_to_char[static_cast<rank_type>(value)];
    }

    //! assign from a character
    constexpr aa27 & assign_char(char_type const c)
    {
        value = char_to_value[c];
        return *this;
    }

    //! return the letter's numeric value or rank in the alphabet
    constexpr rank_type to_rank() const
    {
        return static_cast<rank_type>(value);
    }

    //! assign from a numeric value
    constexpr aa27 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        value = static_cast<internal_type>(c);
        return *this;
    }

    //! The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{27};

    //! @name comparison operators
    //!@{
    constexpr bool operator==(aa27 const & rhs) const
    {
        return value == rhs.value;
    }

    constexpr bool operator!=(aa27 const & rhs) const
    {
        return value != rhs.value;
    }

    constexpr bool operator<(aa27 const & rhs) const
    {
        return value < rhs.value;
    }

    constexpr bool operator>(aa27 const & rhs) const
    {
        return value > rhs.value;
    }

    constexpr bool operator<=(aa27 const & rhs) const
    {
        return value <= rhs.value;
    }

    constexpr bool operator>=(aa27 const & rhs) const
    {
        return value >= rhs.value;
    }
    //!@}

    //! \privatesection
    // conversion tables
    static constexpr char value_to_char[value_size]
    {
        'A',
        'B',
        'C',
        'D',
        'E',
        'F',
        'G',
        'H',
        'I',
        'J',
        'K',
        'L',
        'M',
        'N',
        'O',
        'P',
        'Q',
        'R',
        'S',
        'T',
        'U',
        'V',
        'W',
        'X',
        'Y',
        'Z',
        '*'
    };

    static constexpr internal_type char_to_value[256]
    {
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        //                                              *,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::TERMINATOR, internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        //                      A,                      B,                         C,
        internal_type::UNKNOWN, internal_type::A,       internal_type::B,          internal_type::C,
        //D,                    E,                      F,                         G,
        internal_type::D,       internal_type::E,       internal_type::F,          internal_type::G,
        //H,                    I,                      J,                         K,
        internal_type::H,       internal_type::I,       internal_type::J,          internal_type::K,
        //L,                    M,                      N,                         O,
        internal_type::L,       internal_type::M,       internal_type::N,          internal_type::O,
        //P,                    Q,                      R,                         S,
        internal_type::P,       internal_type::Q,       internal_type::R,          internal_type::S,
        //T,                    U,                      V,                         W,
        internal_type::T,       internal_type::U,       internal_type::V,          internal_type::W,
        //X,                    Y,                      Z
        internal_type::X,       internal_type::Y,       internal_type::Z,          internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        //                      a,                      b,                         c,
        internal_type::UNKNOWN, internal_type::A,       internal_type::B,          internal_type::C,
        //d,                    e,                      f,                         g,
        internal_type::D,       internal_type::E,       internal_type::F,          internal_type::G,
        //h,                    i,                      j,                         k,
        internal_type::H,       internal_type::I,       internal_type::J,          internal_type::K,
        //l,                    m,                      n,                         o,
        internal_type::L,       internal_type::M,       internal_type::N,          internal_type::O,
        //p,                    q,                      r,                         s,
        internal_type::P,       internal_type::Q,       internal_type::R,          internal_type::S,
        //t,                    u,                      v,                         w,
        internal_type::T,       internal_type::U,       internal_type::V,          internal_type::W,
        //x,                    y,                      z
        internal_type::X,       internal_type::Y,       internal_type::Z,          internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN,
        internal_type::UNKNOWN, internal_type::UNKNOWN, internal_type::UNKNOWN,    internal_type::UNKNOWN
    };

};

constexpr aa27 aa27::A{internal_type::A};
constexpr aa27 aa27::B{internal_type::B};
constexpr aa27 aa27::C{internal_type::C};
constexpr aa27 aa27::D{internal_type::D};
constexpr aa27 aa27::E{internal_type::E};
constexpr aa27 aa27::F{internal_type::F};
constexpr aa27 aa27::G{internal_type::G};
constexpr aa27 aa27::H{internal_type::H};
constexpr aa27 aa27::I{internal_type::I};
constexpr aa27 aa27::J{internal_type::J};
constexpr aa27 aa27::K{internal_type::K};
constexpr aa27 aa27::L{internal_type::L};
constexpr aa27 aa27::M{internal_type::M};
constexpr aa27 aa27::N{internal_type::N};
constexpr aa27 aa27::O{internal_type::O};
constexpr aa27 aa27::P{internal_type::P};
constexpr aa27 aa27::Q{internal_type::Q};
constexpr aa27 aa27::R{internal_type::R};
constexpr aa27 aa27::S{internal_type::S};
constexpr aa27 aa27::T{internal_type::T};
constexpr aa27 aa27::U{internal_type::U};
constexpr aa27 aa27::V{internal_type::V};
constexpr aa27 aa27::W{internal_type::W};
constexpr aa27 aa27::X{internal_type::X};
constexpr aa27 aa27::Y{internal_type::Y};
constexpr aa27 aa27::Z{internal_type::Z};
constexpr aa27 aa27::TERMINATOR{internal_type::TERMINATOR};
constexpr aa27 aa27::UNKNOWN{aa27::X};

#ifndef NDEBUG
static_assert(alphabet_concept<aa27>);
#endif

}
