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

/*!\file alphabet/nucleotide/nucl16.hpp
 * \ingroup alphabet
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains seqan3::nucl16, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/concept.hpp>

// ------------------------------------------------------------------
// nucl16
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The 16 letter DNA alphabet, containing all IUPAC smybols.
 * \ingroup alphabet
 *
 * \details
 * Note that in contrast to seqan3::dna4, seqan3::rna4, seqan3::dna5 and seqan3::rna5
 * the letters 'T' and 'U' are distinct values in this alphabet.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     nucl16 my_letter{nucl16::A};
 *     // doesn't work:
 *     // nucl16 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to N internally
 *     if (my_letter.to_char() == 'N')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct nucl16
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
    constexpr nucl16 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr nucl16 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{16};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(nucl16 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(nucl16 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(nucl16 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(nucl16 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(nucl16 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(nucl16 const & rhs) const noexcept
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

    //!\brief Value to char conversion table.
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

    //!\brief Char to value conversion table.
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
public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
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

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::nucl16 is defined as being a nucleotide alphabet.
template <>
struct is_nucleotide<nucl16> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::nucl16>);
static_assert(seqan3::nucleotide_concept<seqan3::nucl16>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{
    /*!\name Alphabet aliases
     * \{
     * \brief Other names (typedefs) for seqan3::nucl16
     * \relates nucl16
     */
    using dna16 = nucl16;
    using rna16 = nucl16;
    using dna = nucl16;
    using rna = nucl16;
    //!\}

} // namespace seqan3


// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::nucl16.
//!\relates nucl16
using nucl16_vector = std::vector<nucl16>;


/*!\brief Alias for an std::basic_string of seqan3::nucl16.
 * \relates nucl16
 *
 * \attention
 * Note that we recommend using seqan3::nucl16_vector instead of nucl16_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using nucl16_string = std::basic_string<nucl16, std::char_traits<nucl16>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief nucl16 literal
 * \relates seqan3::nucl16
 * \returns seqan3::nucl16_vector
 *
 * You can use this string literal to easily assign to nucl16_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // nucl16_vector foo{"ACGTTA"};
 *     // nucl16_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     nucl16_vector foo{"ACGTTA"_nucl16};
 *     nucl16_vector bar = "ACGTTA"_nucl16;
 *     auto bax = "ACGTTA"_nucl16;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline nucl16_vector operator "" _nucl16(const char * s, std::size_t n)
{
    nucl16_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief nucl16 string literal
 * \relates seqan3::nucl16
 * \returns seqan3::nucl16_string
 *
 * You can use this string literal to easily assign to nucl16_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // nucl16_string foo{"ACGTTA"};
 *     // nucl16_string bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     nucl16_string foo{"ACGTTA"_nucl16s};
 *     nucl16_string bar = "ACGTTA"_nucl16s;
 *     auto bax = "ACGTTA"_nucl16s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::nucl16_string and consider using the \link operator""_nucl16 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline nucl16_string operator "" _nucl16s(const char * s, std::size_t n)
{
    nucl16_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal

