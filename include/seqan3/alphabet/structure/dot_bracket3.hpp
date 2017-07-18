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

#include <seqan3/alphabet/structure/concept.hpp>

// ------------------------------------------------------------------
// db3
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The three letter RNA structure alphabet of the characters ".()".
 * \ingroup structure
 *
 * \details
 * The brackets denote RNA base pair interactions. Every left bracket (BL) must have a corresponding right bracket (BR).
 * Pseudoknots cannot be expressed in this format. A dot (.) represents a character that is not paired (NP).
 *
 *~~~~~~~~~~~~~~~{.cpp}
 * GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
 * (((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).
 *~~~~~~~~~~~~~~~
 */

struct db3
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
    static const db3 NP; // not paired
    static const db3 BL; // bracket left
    static const db3 BR; // bracket right
    static const db3 NA; // not available
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
    constexpr db3 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr db3 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{3};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(db3 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(db3 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(db3 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(db3 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(db3 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(db3 const & rhs) const noexcept
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
        UNPAIRED,
        BRACKETL,
        BRACKETR,
        UNKNOWN = UNPAIRED
    };

    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        '.',
        '(',
        ')'
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
            ret['.'] = in_t::UNPAIRED;
            ret['('] = in_t::BRACKETL;
            ret[')'] = in_t::BRACKETR;

            // iupac characters are implicitly "UNKNOWN"
            return ret;
        } ()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr db3 db3::NP{internal_type::UNPAIRED};
constexpr db3 db3::BL{internal_type::BRACKETL};
constexpr db3 db3::BR{internal_type::BRACKETR};
constexpr db3 db3::NA{internal_type::UNKNOWN};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::db3 is defined as being a structure alphabet.
//!\ingroup structure
template <>
struct is_structure<db3> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::db3>);
static_assert(seqan3::structure_concept<seqan3::db3>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::db3.
//!\relates db3
using db3_vector = std::vector<db3>;


/*!\brief Alias for an std::basic_string of seqan3::db3.
 * \relates db3
 *
 * \attention
 * Note that we recommend using seqan3::db3_vector instead of db3_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using db3_string = std::basic_string<db3, std::char_traits<db3>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief db3 literal
 * \relates seqan3::db3
 * \returns seqan3::db3_vector
 *
 * You can use this string literal to easily assign to db3_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     using namespace seqan3::literal;
 *     db3_vector foo{".(..)."_db3};
 *     db3_vector bar = ".(..)."_db3;
 *     auto bax = ".(..)."_db3;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline db3_vector operator "" _db3(const char * s, std::size_t n)
{
    db3_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief db3 string literal
 * \relates seqan3::db3
 * \returns seqan3::db3_string
 *
 * You can use this string literal to easily assign to db3_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     using namespace seqan3::literal;
 *     db3_string foo{".(..)."_db3s};
 *     db3_string bar = ".(..)."_db3s;
 *     auto bax = ".(..)."_db3s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::db3_string and consider using the \link operator""_db3 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline db3_string operator "" _db3s(const char * s, std::size_t n)
{
    db3_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
