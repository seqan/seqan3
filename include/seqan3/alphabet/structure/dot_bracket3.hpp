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
 * \brief Contains the dot bracket format for RNA structure.
 */

#pragma once

#include <cassert>
#include <string>
#include <vector>

#include <seqan3/alphabet/structure/rna_structure_concept.hpp>

// ------------------------------------------------------------------
// dot_bracket3
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The three letter RNA structure alphabet of the characters ".()".
 * \ingroup structure
 *
 * \details
 * The brackets denote RNA base pair interactions. Every left bracket must have a corresponding right bracket.
 * Pseudoknots cannot be expressed in this format. A dot (.) represents a character that is not paired.
 *
 *```console
 *     GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
 *     (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).
 *```
 *
 * \par Usage
 * The following code example creates a dot_bracket3 vector, modifies it, and prints the result to stdout.
 * ```cpp
 *     // create vector
 *     std::vector<dot_bracket3> vec{dot_bracket3::UNPAIRED, dot_bracket3::PAIR_CLOSE, dot_bracket3::PAIR_CLOSE};
 *     // modify and print
 *     vec[1] = dot_bracket3::PAIR_OPEN;
 *     for (dot_bracket3 chr : vec)
 *         std::cout << chr;  // .()
 * ```
 */

struct dot_bracket3
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
    static const dot_bracket3 UNPAIRED;
    static const dot_bracket3 PAIR_OPEN;
    static const dot_bracket3 PAIR_CLOSE;
    static const dot_bracket3 UNKNOWN;
    //!\}

    /*!\name Read functions
     * \{
     *
     * \brief Get the letter as a character of char_type.
     * \returns The character representation of this dot_bracket3 letter.
     */
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    /*!\brief Get the letter's numeric value or rank in the alphabet.
     * \returns The numeric representation of this dot_bracket3 letter.
     */
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }
    //!\}

    /*!\name Write functions
     * \{
     *
     * \brief Assign from a character.
     * \param chr The character that is assigned.
     * \returns The resulting dot_bracket3 character.
     */
    constexpr dot_bracket3 & assign_char(char_type const chr) noexcept
    {
        _value = char_to_value[chr];
        return *this;
    }

    /*!\brief Assign from a numeric value.
     * \param rnk The rank value that is assigned.
     * \returns The resulting dot_bracket3 character.
     */
    constexpr dot_bracket3 & assign_rank(rank_type const rnk)
    {
        assert(rnk < value_size);
        _value = static_cast<internal_type>(rnk);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{3};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(dot_bracket3 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(dot_bracket3 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(dot_bracket3 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(dot_bracket3 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(dot_bracket3 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(dot_bracket3 const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}

    //!\name RNA structure properties
    //!\{
    constexpr bool is_pair_open() const noexcept
    {
        return _value == internal_type::PAIR_OPEN;
    }

    constexpr bool is_pair_close() const noexcept
    {
        return _value == internal_type::PAIR_CLOSE;
    }

    constexpr bool is_unpaired() const noexcept
    {
        return _value == internal_type::UNPAIRED;
    }

    //!\brief The ability of the alphabet to represent pseudoknots, i.e. crossing interactions.
    static constexpr bool pseudoknot_support{false};
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
        UNPAIRED,
        PAIR_OPEN,
        PAIR_CLOSE,
        UNKNOWN = UNPAIRED
    };

    //!\brief Value-to-char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        '.',
        '(',
        ')'
    };

    //!\brief Char-to-value conversion table.
    static constexpr std::array<internal_type, 256> char_to_value
    {
        [] () constexpr
        {
            std::array<internal_type, 256> rank_table{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (internal_type & rnk : rank_table)
                rnk = internal_type::UNKNOWN;

            // canonical
            rank_table['.'] = internal_type::UNPAIRED;
            rank_table['('] = internal_type::PAIR_OPEN;
            rank_table[')'] = internal_type::PAIR_CLOSE;

            return rank_table;
        } ()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr dot_bracket3 dot_bracket3::UNPAIRED{internal_type::UNPAIRED};
constexpr dot_bracket3 dot_bracket3::PAIR_OPEN{internal_type::PAIR_OPEN};
constexpr dot_bracket3 dot_bracket3::PAIR_CLOSE{internal_type::PAIR_CLOSE};
constexpr dot_bracket3 dot_bracket3::UNKNOWN{internal_type::UNKNOWN};

} // namespace seqan3

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::dot_bracket3>);
static_assert(seqan3::rna_structure_concept<seqan3::dot_bracket3>);
#endif

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief dot_bracket3 literal
 * \relates seqan3::dot_bracket3
 * \returns std::vector<seqan3::dot_bracket3>
 *
 * You can use this string literal to easily assign to a vector of dot_bracket3 characters:
 *
 *```.cpp
 *     using namespace seqan3::literal;
 *     std::vector<dot_bracket3> foo{".(..)."_db3};
 *     std::vector<dot_bracket3> bar = ".(..)."_db3;
 *     auto bax = ".(..)."_db3;
 *```
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */
inline std::vector<dot_bracket3> operator "" _db3(const char * str, std::size_t len)
{
    std::vector<dot_bracket3> vec;
    vec.resize(len);

    for (size_t idx = 0u; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

/*!\brief dot_bracket3 string literal
 * \relates seqan3::dot_bracket3
 * \returns std::basic_string<seqan3::dot_bracket3, std::char_traits<seqan3::dot_bracket3>>
 *
 * You can use this string literal to easily assign to a string of dot_bracket3 characters:
 *
 *```.cpp
 *     using namespace seqan3::literal;
 *     using string_t = std::basic_string<dot_bracket3, std::char_traits<dot_bracket3>>;
 *     string_t foo{".(..)."_db3s};
 *     string_t bar = ".(..)."_db3s;
 *     auto bax = ".(..)."_db3s;
 *```
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */
inline std::basic_string<dot_bracket3, std::char_traits<dot_bracket3>> operator "" _db3s(const char * str,
                                                                                                  std::size_t len)
{
    std::basic_string<dot_bracket3, std::char_traits<dot_bracket3>> db3str;
    db3str.resize(len);

    for (size_t idx = 0u; idx < len; ++idx)
        db3str[idx].assign_char(str[idx]);

    return db3str;
}

} // namespace seqan3::literal
