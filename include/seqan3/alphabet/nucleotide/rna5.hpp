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

/*!\file alphabet/nucleotide/rna5.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::rna5, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>

// ------------------------------------------------------------------
// rna5
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The five letter RNA alphabet of A,C,G,U and the unknown character N.
 * \ingroup alphabet
 *
 * \details
 * This alphabet inherits from seqan3::dna5 and is guaranteed to have the same internal representation of
 * data. The only difference is that it prints 'U' on character conversion instead of 'T'. You assign
 * between values of seqan3::dna5 and seqan3::rna5.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     rna5 my_letter{rna5::A};
 *     // doesn't work:
 *     // rna5 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to A internally
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct rna5 : public dna5
{
    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.*
     */
    //!\{
    static const rna5 A;
    static const rna5 C;
    static const rna5 G;
    static const rna5 T;
    static const rna5 N;
    static const rna5 U;
    static const rna5 UNKNOWN;
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type.
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    constexpr rna5 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr rna5 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }

    //!\brief Assign from seqan3::dna5.
    constexpr rna5 & operator=(dna5 const in) noexcept
    {
        _value = in._value;
        return *this;
    }
    //!\}

private:
    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'U',
        'N'
    };
};

constexpr rna5 rna5::A{internal_type::A};
constexpr rna5 rna5::C{internal_type::C};
constexpr rna5 rna5::G{internal_type::G};
constexpr rna5 rna5::U{internal_type::U};
constexpr rna5 rna5::N{internal_type::N};
constexpr rna5 rna5::T{rna5::U};
constexpr rna5 rna5::UNKNOWN{rna5::N};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::rna5 is defined as being a nucleotide alphabet.
template <>
struct is_nucleotide<rna5> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::rna5>);
static_assert(seqan3::nucleotide_concept<seqan3::rna5>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::rna5.
//!\relates rna5
using rna5_vector = std::vector<rna5>;


/*!\brief Alias for an std::basic_string of seqan3::rna5.
 * \relates rna5
 *
 * \attention
 * Note that we recommend using seqan3::rna5_vector instead of rna5_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using rna5_string = std::basic_string<rna5, std::char_traits<rna5>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief rna5 literal
 * \relates seqan3::rna5
 * \returns seqan3::rna5_vector
 *
 * You can use this string literal to easily assign to rna5_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna5_vector foo{"ACGTTA"};
 *     // rna5_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     rna5_vector foo{"ACGTTA"_rna5};
 *     rna5_vector bar = "ACGTTA"_rna5;
 *     auto bax = "ACGTTA"_rna5;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline rna5_vector operator "" _rna5(const char * s, std::size_t n)
{
    rna5_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief rna5 string literal
 * \relates seqan3::rna5
 * \returns seqan3::rna5_string
 *
 * You can use this string literal to easily assign to rna5_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna5_string foo{"ACGTTA"};
 *     // rna5_string bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     rna5_string foo{"ACGTTA"_rna5s};
 *     rna5_string bar = "ACGTTA"_rna5s;
 *     auto bax = "ACGTTA"_rna5s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::rna5_string and consider using the \link operator""_rna5 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline rna5_string operator "" _rna5s(const char * s, std::size_t n)
{
    rna5_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal

