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

/*!\file alphabet/nucleotide/dna4.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::dna4, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/concept.hpp>

// ------------------------------------------------------------------
// dna4
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The four letter DNA alphabet of A,C,G,T.
 * \ingroup alphabet
 *
 * \details
 * Note that you can assign 'U' as a character to dna4 and it will silently
 * be converted to 'T'.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     dna4 my_letter{dna4::A};
 *     // doesn't work:
 *     // dna4 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to A internally
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct dna4
{
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = uint8_t;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface. *Don't worry about the `internal_type`.*
     */
    //!\{
    static const dna4 A;
    static const dna4 C;
    static const dna4 G;
    static const dna4 T;
    static const dna4 U;
    static const dna4 UNKNOWN;
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
    constexpr dna4 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr dna4 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{4};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(dna4 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(dna4 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(dna4 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(dna4 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(dna4 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(dna4 const & rhs) const noexcept
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
        C,
        G,
        T,
        U = T,
        UNKNOWN = A
    };

    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'T'
    };

    //!\brief Char to value conversion table.
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
public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr dna4 dna4::A{internal_type::A};
constexpr dna4 dna4::C{internal_type::C};
constexpr dna4 dna4::G{internal_type::G};
constexpr dna4 dna4::T{internal_type::T};
constexpr dna4 dna4::U{dna4::T};
constexpr dna4 dna4::UNKNOWN{dna4::A};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::dna4 is defined as being a nucleotide alphabet.
template <>
struct is_nucleotide<dna4> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::dna4>);
static_assert(seqan3::nucleotide_concept<seqan3::dna4>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::dna4.
//!\relates dna4
using dna4_vector = std::vector<dna4>;


/*!\brief Alias for an std::basic_string of seqan3::dna4.
 * \relates dna4
 *
 * \attention
 * Note that we recommend using seqan3::dna4_vector instead of dna4_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using dna4_string = std::basic_string<dna4, std::char_traits<dna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief dna4 literal
 * \relates seqan3::dna4
 * \returns seqan3::dna4_vector
 *
 * You can use this string literal to easily assign to dna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna4_vector foo{"ACGTTA"};
 *     // dna4_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     dna4_vector foo{"ACGTTA"_dna4};
 *     dna4_vector bar = "ACGTTA"_dna4;
 *     auto bax = "ACGTTA"_dna4;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All user-defined literals are in the namespace seqan3::literal!
 */

inline dna4_vector operator "" _dna4(const char * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief dna4 string literal
 * \relates seqan3::dna4
 * \returns seqan3::dna4_string
 *
 * You can use this string literal to easily assign to dna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna4_string foo{"ACGTTA"};
 *     // dna4_string bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     dna4_string foo{"ACGTTA"_dna4s};
 *     dna4_string bar = "ACGTTA"_dna4s;
 *     auto bax = "ACGTTA"_dna4s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::dna4_string and consider using the \link operator""_dna4 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline dna4_string operator "" _dna4s(const char * s, std::size_t n)
{
    dna4_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
