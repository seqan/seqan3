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
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains the dssp format for protein structure.
 */

#pragma once

#include <cassert>
#include <vector>

#include <seqan3/alphabet/concept.hpp>

// ------------------------------------------------------------------
// dssp9
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The protein structure alphabet of the characters "HGIEBTSCX".
 * \ingroup structure
 *
 * \details
 * The DSSP annotation links structure elements to protein sequences.
 * Originally created with 7 letters as a file format for the DSSP program (http://www.cmbi.ru.nl/dssp.html),
 * it is also used in the stockholm file format for structure alignments, extended by the characters C and X
 * (https://en.wikipedia.org/wiki/Stockholm_format).
 *
 * The letter abbreviations are as follows:
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
 *
 * \par Usage
 * The following code example creates a dssp9 vector, modifies it, and prints the result to stdout.
 * ```cpp
 *     // create vector
 *     std::vector<dssp9> vec{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
 *     // modify and print
 *     vec[1] = dssp9::C;
 *     for (dssp9 chr : vec)
 *         std::cout << chr;  // ECHHTG
 * ```
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

    //!\name Read functions
    //!\{

    /*!\brief Get the letter as a character of char_type.
     * \returns The character representation of this dssp9 letter.
     */
    constexpr char_type to_char() const noexcept { return value_to_char[static_cast<rank_type>(_value)]; }

    /*!\brief Get the letter's numeric value or rank in the alphabet.
     * \returns The numeric representation of this dssp9 letter.
     */
    constexpr rank_type to_rank() const noexcept { return static_cast<rank_type>(_value); }
    //!\}

    //!\name Write functions
    //!\{

    /*!\brief Assign from a character.
     * \param chr The character that is assigned.
     * \returns The resulting dssp9 character.
     */
    constexpr dssp9 & assign_char(char_type const chr) noexcept
    {
        _value = char_to_value[chr];
        return *this;
    }

    /*!\brief Assign from a numeric value.
     * \param rnk The rank value that is assigned.
     * \returns The resulting dssp9 character.
     */
    constexpr dssp9 & assign_rank(rank_type const rnk)
    {
        assert(rnk < value_size);
        _value = static_cast<internal_type>(rnk);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{ 9 };

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(dssp9 const & rhs) const noexcept { return _value == rhs._value; }

    constexpr bool operator!=(dssp9 const & rhs) const noexcept { return _value != rhs._value; }

    constexpr bool operator<(dssp9 const & rhs) const noexcept { return _value < rhs._value; }

    constexpr bool operator>(dssp9 const & rhs) const noexcept { return _value > rhs._value; }

    constexpr bool operator<=(dssp9 const & rhs) const noexcept { return _value <= rhs._value; }

    constexpr bool operator>=(dssp9 const & rhs) const noexcept { return _value >= rhs._value; }
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
        H,
        B,
        E,
        G,
        I,
        T,
        S,
        C,
        X
    };

    //!\brief Value-to-char conversion table.
    static constexpr char_type value_to_char[value_size]{ 'H', 'B', 'E', 'G', 'I', 'T', 'S', 'C', 'X' };

    //!\brief Char-to-value conversion table.
    static constexpr std::array<internal_type, 256> char_to_value{
        []() constexpr { std::array<internal_type, 256> rank_table{};

    // initialize with X (std::array::fill unfortunately not constexpr)
    for (internal_type & rnk : rank_table) rnk = internal_type::X;

    // canonical
    rank_table['H'] = internal_type::H;
    rank_table['B'] = internal_type::B;
    rank_table['E'] = internal_type::E;
    rank_table['G'] = internal_type::G;
    rank_table['I'] = internal_type::I;
    rank_table['T'] = internal_type::T;
    rank_table['S'] = internal_type::S;
    rank_table['C'] = internal_type::C;
    rank_table['X'] = internal_type::X;

    return rank_table;
}()
};

public:
//!\privatesection
//!\brief The data member.
internal_type _value;
//!\publicsection
}
;

constexpr dssp9 dssp9::H{ internal_type::H };
constexpr dssp9 dssp9::B{ internal_type::B };
constexpr dssp9 dssp9::E{ internal_type::E };
constexpr dssp9 dssp9::G{ internal_type::G };
constexpr dssp9 dssp9::I{ internal_type::I };
constexpr dssp9 dssp9::T{ internal_type::T };
constexpr dssp9 dssp9::S{ internal_type::S };
constexpr dssp9 dssp9::C{ internal_type::C };
constexpr dssp9 dssp9::X{ internal_type::X };

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief dssp9 literal
 * \relates seqan3::dssp9
 * \returns std::vector<seqan3::dssp9>
 *
 * You can use this string literal to easily assign to a vector of dssp9 characters:
 *
 *```.cpp
 *     using namespace seqan3::literal;
 *     std::vector<dssp9> foo{"EHHHHT"_dssp9};
 *     std::vector<dssp9> bar = "EHHHHT"_dssp9;
 *     auto bax = "EHHHHT"_dssp9;
 *```
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */
inline std::vector<dssp9> operator""_dssp9(const char * str, std::size_t len)
{
    std::vector<dssp9> vec;
    vec.resize(len);

    for (size_t idx = 0u; idx < len; ++idx) vec[idx].assign_char(str[idx]);

    return vec;
}

} // namespace seqan3::literal
