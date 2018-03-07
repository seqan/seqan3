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
 * \ingroup nucleotide
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::rna15, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>

// ------------------------------------------------------------------
// rna15
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The 15 letter RNA alphabet, containing all IUPAC smybols minus the gap.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_concept
 *
 * \details
 * This alphabet inherits from seqan3::dna15 and is guaranteed to have the same internal representation of
 * data. The only difference is that it prints 'U' on character conversion instead of 'T'. You can implicitly
 * convert between values of seqan3::dna15 and seqan3::rna15.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     rna15 my_letter{rna15::A};
 *     // doesn't work:
 *     // rna15 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to N internally
 *     if (my_letter.to_char() == 'N')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct rna15 : public dna15
{
    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.*
     */
    //!\{
    static const rna15 A;
    static const rna15 B;
    static const rna15 C;
    static const rna15 D;
    static const rna15 G;
    static const rna15 H;
    static const rna15 K;
    static const rna15 M;
    static const rna15 N;
    static const rna15 R;
    static const rna15 S;
    static const rna15 T;
    static const rna15 U;
    static const rna15 V;
    static const rna15 W;
    static const rna15 Y;
    static const rna15 UNKNOWN;
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\copydoc seqan3::dna4::to_char
    constexpr char_type to_char() const noexcept { return value_to_char[static_cast<rank_type>(_value)]; }

    //!\copydoc seqan3::dna4::complement
    constexpr rna15 complement() const noexcept { return dna15::complement(); }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\copydoc seqan3::dna4::assign_char
    constexpr rna15 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\copydoc seqan3::dna4::assign_rank
    constexpr rna15 & assign_rank(rank_type const c)
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
        return other_nucl_type{ _value };
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

  protected:
    //!\privatesection
    //!\copydoc seqan3::dna4::value_to_char
    static constexpr char_type value_to_char[value_size]{ 'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M',
                                                          'N', 'R', 'S', 'U', 'V', 'W', 'Y' };
};

constexpr rna15 rna15::A{ internal_type::A };
constexpr rna15 rna15::B{ internal_type::B };
constexpr rna15 rna15::C{ internal_type::C };
constexpr rna15 rna15::D{ internal_type::D };
constexpr rna15 rna15::G{ internal_type::G };
constexpr rna15 rna15::H{ internal_type::H };
constexpr rna15 rna15::K{ internal_type::K };
constexpr rna15 rna15::M{ internal_type::M };
constexpr rna15 rna15::N{ internal_type::N };
constexpr rna15 rna15::R{ internal_type::R };
constexpr rna15 rna15::S{ internal_type::S };
constexpr rna15 rna15::U{ internal_type::U };
constexpr rna15 rna15::T{ rna15::U };
constexpr rna15 rna15::V{ internal_type::V };
constexpr rna15 rna15::W{ internal_type::W };
constexpr rna15 rna15::Y{ internal_type::Y };
constexpr rna15 rna15::UNKNOWN{ rna15::N };

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::rna15.
//!\relates rna15
using rna15_vector = std::vector<rna15>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief rna15 literal
 * \relates seqan3::rna15
 * \returns seqan3::rna15_vector
 *
 * You can use this string literal to easily assign to rna15_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // rna15_vector foo{"ACGTTA"};
 *     // rna15_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     rna15_vector foo{"ACGTTA"_rna15};
 *     rna15_vector bar = "ACGTTA"_rna15;
 *     auto bax = "ACGTTA"_rna15;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline rna15_vector operator""_rna15(const char * s, std::size_t n)
{
    rna15_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i) r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
