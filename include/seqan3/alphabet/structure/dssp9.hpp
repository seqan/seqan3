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

/*!\file
 * \ingroup structure
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains the dssp format for protein structure.
 */

#pragma once

#include <cassert>
#include <string>
#include <vector>

#include <seqan3/alphabet/structure/concept.hpp>

// ------------------------------------------------------------------
// dssp9
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The protein structure alphabet of the characters "HGIEBTSCX".
 * \ingroup structure
 *
 * \details
 * Format description regarding to the DSSP website (\link http://www.cmbi.ru.nl/dssp.html \endlink)
 * and the Stockholm format description (\link https://en.wikipedia.org/wiki/Stockholm_format \endlink)
 *
 * H = alpha helix
 * B = beta bridge
 * E = strand
 * G = helix-3
 * I = helix-5
 * T = turn
 * S = bend
 * C = coil/loop
 * X = unknown
 */

struct dssp9
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
    static const dssp9 H;
    static const dssp9 B;
    static const dssp9 E;
    static const dssp9 G;
    static const dssp9 I;
    static const dssp9 T;
    static const dssp9 S;
    static const dssp9 C;
    static const dssp9 X;
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
    constexpr dssp9 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr dssp9 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{9};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(dssp9 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(dssp9 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(dssp9 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(dssp9 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(dssp9 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(dssp9 const & rhs) const noexcept
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
        H, B, E, G, I, T, S, C, X
    };

    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        'H', 'B', 'E', 'G', 'I', 'T', 'S', 'C', 'X'
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
            c = in_t::X;

        // canonical
        ret['H'] = in_t::H;
        ret['B'] = in_t::B;
        ret['E'] = in_t::E;
        ret['G'] = in_t::G;
        ret['I'] = in_t::I;
        ret['T'] = in_t::T;
        ret['S'] = in_t::S;
        ret['C'] = in_t::C;
        ret['X'] = in_t::X;

        return ret;
    } ()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr dssp9 dssp9::H{internal_type::H};
constexpr dssp9 dssp9::B{internal_type::B};
constexpr dssp9 dssp9::E{internal_type::E};
constexpr dssp9 dssp9::G{internal_type::G};
constexpr dssp9 dssp9::I{internal_type::I};
constexpr dssp9 dssp9::T{internal_type::T};
constexpr dssp9 dssp9::S{internal_type::S};
constexpr dssp9 dssp9::C{internal_type::C};
constexpr dssp9 dssp9::X{internal_type::X};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::dssp9 is defined as being a structure alphabet.
//!\ingroup structure
template <>
struct is_structure<dssp9> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::dssp9>);
static_assert(seqan3::structure_concept<seqan3::dssp9>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::dssp9.
//!\relates dssp9
using dssp9_vector = std::vector<dssp9>;


/*!\brief Alias for an std::basic_string of seqan3::dssp9.
 * \relates dssp9
 *
 * \attention
 * Note that we recommend using seqan3::dssp9_vector instead of dssp9_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using dssp9_string = std::basic_string<dssp9, std::char_traits<dssp9>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief dssp9 literal
 * \relates seqan3::dssp9
 * \returns seqan3::dssp9_vector
 *
 * You can use this string literal to easily assign to dssp9_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     using namespace seqan3::literal;
 *     dssp9_vector foo{"EHHHHT"_dssp9};
 *     dssp9_vector bar = "EHHHHT"_dssp9;
 *     auto bax = "EHHHHT"_dssp9;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline dssp9_vector operator "" _dssp9(const char * s, std::size_t n)
{
    dssp9_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief dssp9 string literal
 * \relates seqan3::dssp9
 * \returns seqan3::dssp9_string
 *
 * You can use this string literal to easily assign to dssp9_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     using namespace seqan3::literal;
 *     dssp9_string foo{"EHHHHT"_dssp9s};
 *     dssp9_string bar = "EHHHHT"_dssp9s;
 *     auto bax = "EHHHHT"_dssp9s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::dssp9_string and consider using the \link operator""_dssp9 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline dssp9_string operator "" _dssp9s(const char * s, std::size_t n)
{
    dssp9_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
