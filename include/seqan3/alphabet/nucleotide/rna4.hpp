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

/*!\file alphabet/nucleotide/rna4.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::rna4, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

// ------------------------------------------------------------------
// rna4
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The four letter RNA alphabet of A,C,G,U.
 * \ingroup alphabet
 *
 * \details
 * This alphabet inherits from seqan3::dna4 and is guaranteed to have the same internal representation of
 * data. The only difference is that it prints 'U' on character conversion instead of 'T'. You assign
 * between values of seqan3::dna4 and seqan3::rna4.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     rna4 my_letter{rna4::A};
 *     // doesn't work:
 *     // rna4 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to A internally
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct rna4 : public dna4
{
    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.*
     */
    //!\{
    static const rna4 A;
    static const rna4 C;
    static const rna4 G;
    static const rna4 U;
    static const rna4 T;
    static const rna4 UNKNOWN;
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
    constexpr rna4 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr rna4 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }

    //!\brief Assign from seqan3::dna4.
    constexpr rna4 & operator=(dna4 const in) noexcept
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
        'U'
    };
};

constexpr rna4 rna4::A{internal_type::A};
constexpr rna4 rna4::C{internal_type::C};
constexpr rna4 rna4::G{internal_type::G};
constexpr rna4 rna4::U{internal_type::U};
constexpr rna4 rna4::T{rna4::U};
constexpr rna4 rna4::UNKNOWN{rna4::A};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::rna4 is defined as being a nucleotide alphabet.
template <>
struct is_nucleotide<rna4> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::rna4>);
static_assert(seqan3::nucleotide_concept<seqan3::rna4>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::rna4.
//!\relates rna4
using rna4_vector = std::vector<rna4>;


/*!\brief Alias for an std::basic_string of seqan3::rna4.
 * \relates rna4
 *
 * \attention
 * Note that we recommend using seqan3::rna4_vector instead of rna4_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using rna4_string = std::basic_string<rna4, std::char_traits<rna4>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief rna4 literal
 * \relates seqan3::rna4
 * \returns seqan3::rna4_vector
 *
 * You can use this string literal to easily assign to rna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna4_vector foo{"ACGTTA"};
 *     // rna4_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     rna4_vector foo{"ACGTTA"_rna4};
 *     rna4_vector bar = "ACGTTA"_rna4;
 *     auto bax = "ACGTTA"_rna4;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline rna4_vector operator "" _rna4(const char * s, std::size_t n)
{
    rna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief rna4 string literal
 * \relates seqan3::rna4
 * \returns seqan3::rna4_string
 *
 * You can use this string literal to easily assign to rna4_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna4_string foo{"ACGTTA"};
 *     // rna4_string bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     rna4_string foo{"ACGTTA"_rna4s};
 *     rna4_string bar = "ACGTTA"_rna4s;
 *     auto bax = "ACGTTA"_rna4s;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::rna4_string and consider using the \link operator""_rna4 \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline rna4_string operator "" _rna4s(const char * s, std::size_t n)
{
    rna4_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal

