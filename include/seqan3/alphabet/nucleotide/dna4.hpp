// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::dna4, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <type_traits>
#include <vector>

#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

// ------------------------------------------------------------------
// dna4
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The four letter DNA alphabet of A,C,G,T.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_concept
 *
 * \details
 * Note that you can assign 'U' as a character to dna4 and it will silently
 * be converted to 'T'.
 *
 * \snippet test/snippet/alphabet/nucleotide/dna4.cpp code
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
    /*!\brief Return the letter as a character of char_type.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::to_char() requirement via the seqan3::to_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    /*!\brief Return the letter's numeric value or rank in the alphabet.
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::to_rank() requirement via the seqan3::to_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }

    /*!\brief Return the complement of the letter.
     *
     * \details
     *
     * See \ref nucleotide for the actual values.
     *
     * Satisfies the seqan3::nucleotide_concept::complement() requirement via the seqan3::complement() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr dna4 complement() const noexcept
    {
        return complement_table[to_rank()];
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a character.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::assign_char() requirement via the seqan3::assign_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr dna4 & assign_char(char_type const c) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(c)];
        return *this;
    }

    /*!\brief Assign from a numeric value.
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::assign_rank() requirement via the seqan3::assign_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr dna4 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{4};

    /*!\name Conversion operators
     * \{
     */
    //!\brief Implicit conversion between dna* and rna* of the same size.
    //!\tparam other_nucl_type The type to convert to; must satisfy seqan3::nucleotide_concept and have the same \link value_size \endlink.
    template <typename other_nucl_type>
    //!\cond
        requires nucleotide_concept<other_nucl_type> && value_size == alphabet_size_v<other_nucl_type>
    //!\endcond
    constexpr operator other_nucl_type() const noexcept
    {
        return other_nucl_type{_value};
    }

    //!\brief Explicit conversion to any other nucleotide alphabet (via char representation).
    //!\tparam other_nucl_type The type to convert to; must satisfy seqan3::nucleotide_concept.
    template <typename other_nucl_type>
    //!\cond
        requires nucleotide_concept<other_nucl_type>
    //!\endcond
    explicit constexpr operator other_nucl_type() const noexcept
    {
        return detail::convert_through_char_representation<other_nucl_type, std::decay_t<decltype(*this)>>[to_rank()];
    }
    //!\}

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
    static constexpr std::array<internal_type, 256> char_to_value
    {
        [] () constexpr
        {
            using in_t = internal_type;
            std::array<in_t, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = in_t::UNKNOWN;

            // canonical
            ret['A'] = in_t::A; ret['a'] = in_t::A;
            ret['C'] = in_t::C; ret['c'] = in_t::C;
            ret['G'] = in_t::G; ret['g'] = in_t::G;
            ret['T'] = in_t::T; ret['t'] = in_t::T;
            ret['U'] = in_t::U; ret['u'] = in_t::U;

            // iupac characters get special treatment, because there is no N
            ret['R'] = in_t::A; ret['r'] = in_t::A; // or G
            ret['Y'] = in_t::C; ret['y'] = in_t::C; // or T
            ret['S'] = in_t::C; ret['s'] = in_t::C; // or G
            ret['W'] = in_t::A; ret['w'] = in_t::A; // or T
            ret['K'] = in_t::G; ret['k'] = in_t::G; // or T
            ret['M'] = in_t::A; ret['m'] = in_t::A; // or T
            ret['B'] = in_t::C; ret['b'] = in_t::C; // or G or T
            ret['D'] = in_t::A; ret['d'] = in_t::A; // or G or T
            ret['H'] = in_t::A; ret['h'] = in_t::A; // or C or T
            ret['V'] = in_t::A; ret['v'] = in_t::A; // or C or G

            return ret;
        }()
    };

    //!\brief The complement table.
    static const std::array<dna4, value_size> complement_table;

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

constexpr std::array<dna4, dna4::value_size> dna4::complement_table
{
    dna4::T,    // complement of dna4::A
    dna4::G,    // complement of dna4::C
    dna4::C,    // complement of dna4::G
    dna4::A     // complement of dna4::T
};

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::dna4.
//!\relates dna4
using dna4_vector = std::vector<dna4>;

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
 * \snippet test/snippet/alphabet/nucleotide/dna4.cpp operator""_dna4
 *
 * \attention
 * All user-defined literals are in the namespace seqan3::literal!
 */

inline dna4_vector operator""_dna4(const char * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
