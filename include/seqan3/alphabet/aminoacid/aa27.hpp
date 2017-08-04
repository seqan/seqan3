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

/*!\file aa27.hpp
 * \ingroup aminoacid
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains seqan3::aa27, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\brief The twenty-seven letter amino acid alphabet
 * \ingroup aminoacid
 * \implements seqan3::alphabet_concept
 *
 * \details
 * The alphabet consists of letters A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X,
 * Y, Z, *
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 *```cpp
 *     aa27 my_letter{aa27::A};
 *     // doesn't work:
 *     // aa27 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('?'); // converted to X internally
 *     if (my_letter.to_char() == 'X')
 *        std::cout << "yeah\n"; // "yeah";
 *```
 */

struct aa27
{
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;

    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = uint8_t;


    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.*
     */
    //!\{
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
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type.
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    //!\brief Return the letter's numeric value or rank in the alphabet.
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    constexpr aa27 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr aa27 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{27};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(aa27 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(aa27 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(aa27 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(aa27 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(aa27 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(aa27 const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}

    protected:
    //!\privatesection
    /*!\brief The internal type is a strictly typed enum.
     *
     * This is done to prevent aggregate initialization from numbers and/or chars.
     * It is has the drawback that it also introduces a scope which in turn makes
     * the static "letter values " members necessary.
     */
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

    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
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

    //!\brief Char to value conversion table.
    static constexpr std::array<internal_type, 256> char_to_value
    {
        [] () constexpr
        {
            using in_t = internal_type;
            std::array<in_t, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = in_t::UNKNOWN;

            ret['A'] = in_t::A; ret['a'] = in_t::A;
            ret['B'] = in_t::B; ret['b'] = in_t::B;
            ret['C'] = in_t::C; ret['c'] = in_t::C;
            ret['D'] = in_t::D; ret['d'] = in_t::D;
            ret['E'] = in_t::E; ret['e'] = in_t::E;
            ret['F'] = in_t::F; ret['f'] = in_t::F;
            ret['G'] = in_t::G; ret['g'] = in_t::G;
            ret['H'] = in_t::H; ret['h'] = in_t::H;
            ret['I'] = in_t::I; ret['i'] = in_t::I;
            ret['J'] = in_t::J; ret['j'] = in_t::J;
            ret['K'] = in_t::K; ret['k'] = in_t::K;
            ret['L'] = in_t::L; ret['l'] = in_t::L;
            ret['M'] = in_t::M; ret['m'] = in_t::M;
            ret['N'] = in_t::N; ret['n'] = in_t::N;
            ret['O'] = in_t::O; ret['o'] = in_t::O;
            ret['P'] = in_t::P; ret['p'] = in_t::P;
            ret['Q'] = in_t::Q; ret['q'] = in_t::Q;
            ret['R'] = in_t::R; ret['r'] = in_t::R;
            ret['S'] = in_t::S; ret['s'] = in_t::S;
            ret['T'] = in_t::T; ret['t'] = in_t::T;
            ret['U'] = in_t::U; ret['u'] = in_t::U;
            ret['V'] = in_t::V; ret['v'] = in_t::V;
            ret['W'] = in_t::W; ret['w'] = in_t::W;
            ret['X'] = in_t::X; ret['x'] = in_t::X;
            ret['Y'] = in_t::Y; ret['y'] = in_t::Y;
            ret['Z'] = in_t::Z; ret['z'] = in_t::Z;
            ret['*'] = in_t::TERMINATOR;
            return ret;
        }()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
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

} // namespace seqan3

#ifndef NDEBUGs

#include <seqan3/alphabet/concept.hpp>

static_assert(seqan3::alphabet_concept<seqan3::aa27>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{
//!\brief Alias for an std::vector of seqan3::aa27.
//!\relates aa27
using aa27_vector = std::vector<aa27>;

/*!\brief Alias for an std::basic_string of seqan3::aa27.
 * \relates aa27
 *
 * \attention
 * Note that we recommend using seqan3::aa27_vector instead of aa27_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using aa27_string = std::basic_string<aa27, std::char_traits<aa27>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief aa27 literal
 * \relates seqan3::aa27
 * \returns seqan3::aa27_vector
 *
 * You can use this string literal to easily assign to aa27_vector:
 *
 *```cpp
 *     // these don't work:
 *     // aa27_vector foo{"ABFUYR"};
 *     // aa27_vector bar = "ABFUYR";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     aa27_vector foo{"ABFUYR"_aa27};
 *     aa27_vector bar = "ABFUYR"_aa27;
 *     auto bax = "ABFUYR"_aa27;
 *```
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline aa27_vector operator "" _aa27(const char * s, std::size_t n)
{
    aa27_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief aa27 string literal
 * \relates seqan3::aa27
 * \returns seqan3::aa27_string
 *
 * You can use this string literal to easily assign to aa27_vector:
 *
 *```cpp
 *     // these don't work:
 *     // aa27_string foo{"ABFUYR"};
 *     // aa27_string bar = "ABFUYR";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     aa27_string foo{"ABFUYR"_aa27s};
 *     aa27_string bar = "ABFUYR"_aa27s;
 *     auto bax = "ABFUYR"_aa27s;
 *```
 *
 * Please note the limitations of seqan3::aa27_string and consider using the \link operator""_aa27 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline aa27_string operator "" _aa27s(const char * s, std::size_t n)
{
    aa27_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
