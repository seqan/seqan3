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
 * \brief Contains seqan3::rna4, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <vector>

#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

// ------------------------------------------------------------------
// rna4
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The four letter RNA alphabet of A,C,G,U.
 * \ingroup nucleotide
 *
 * \details
 * This alphabet inherits from seqan3::dna4 and is guaranteed to have the same internal representation of
 * data. The only difference is that it prints 'U' on character conversion instead of 'T'. You assign
 * between values of seqan3::dna4 and seqan3::rna4.
 *
 *\snippet test/snippet/alphabet/nucleotide/rna4.cpp code
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
    //!\copydoc seqan3::dna4::to_char
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    //!\copydoc seqan3::dna4::complement
    constexpr rna4 complement() const noexcept
    {
        return rna4{dna4::complement()};
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\copydoc seqan3::dna4::assign_char
    constexpr rna4 & assign_char(char_type const c) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(c)];
        return *this;
    }

    //!\copydoc seqan3::dna4::assign_rank
    constexpr rna4 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

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

private:
    //!\copydoc seqan3::dna4::value_to_char
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

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::rna4.
//!\relates rna4
using rna4_vector = std::vector<rna4>;

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
 *\snippet test/snippet/alphabet/nucleotide/rna4.cpp operator""_rna4
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline rna4_vector operator""_rna4(const char * s, std::size_t n)
{
    rna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
