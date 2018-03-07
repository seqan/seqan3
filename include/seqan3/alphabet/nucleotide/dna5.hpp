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
 * \author David Heller <david.heller AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::dna5, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <vector>

#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

// ------------------------------------------------------------------
// dna5
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The five letter DNA alphabet of A,C,G,T and the unknown character N.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_concept
 *
 * \details
 * Note that you can assign 'U' as a character to dna5 and it will silently
 * be converted to 'T'.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     dna5 my_letter{dna5::A};
 *     // doesn't work:
 *     // dna5 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to N internally
 *     if (my_letter.to_char() == 'N')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct dna5
{
    //!\copydoc seqan3::dna4::char_type
    using char_type = char;
    //!\copydoc seqan3::dna4::rank_type
    using rank_type = uint8_t;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.*
     */
    //!\{
    static const dna5 A;
    static const dna5 C;
    static const dna5 G;
    static const dna5 T;
    static const dna5 N;
    static const dna5 U;
    static const dna5 UNKNOWN;
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\copydoc seqan3::dna4::to_char
    constexpr char_type to_char() const noexcept { return value_to_char[static_cast<rank_type>(_value)]; }

    //!\copydoc seqan3::dna4::to_rank
    constexpr rank_type to_rank() const noexcept { return static_cast<rank_type>(_value); }

    //!\copydoc seqan3::dna4::complement
    constexpr dna5 complement() const noexcept { return complement_table[to_rank()]; }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\copydoc seqan3::dna4::assign_char
    constexpr dna5 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\copydoc seqan3::dna4::assign_rank
    constexpr dna5 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\copydoc seqan3::dna4::value_size
    static constexpr rank_type value_size{ 5 };

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

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(dna5 const & rhs) const noexcept { return _value == rhs._value; }

    constexpr bool operator!=(dna5 const & rhs) const noexcept { return _value != rhs._value; }

    constexpr bool operator<(dna5 const & rhs) const noexcept { return _value < rhs._value; }

    constexpr bool operator>(dna5 const & rhs) const noexcept { return _value > rhs._value; }

    constexpr bool operator<=(dna5 const & rhs) const noexcept { return _value <= rhs._value; }

    constexpr bool operator>=(dna5 const & rhs) const noexcept { return _value >= rhs._value; }
    //!\}

  protected:
    //!\privatesection
    //!\copydoc seqan3::dna4::internal_type
    enum struct internal_type : rank_type
    {
        A,
        C,
        G,
        T,
        N,
        U = T,
        UNKNOWN = N
    };

    //!\copydoc seqan3::dna4::value_to_char
    static constexpr char_type value_to_char[value_size]{ 'A', 'C', 'G', 'T', 'N' };

    //!\copydoc seqan3::dna4::char_to_value
    static constexpr std::array<internal_type, 256> char_to_value{ []() constexpr { using in_t = internal_type;
    std::array<in_t, 256> ret{};

    // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
    for (auto & c : ret) c = in_t::UNKNOWN;

    // canonical
    ret['A'] = in_t::A;
    ret['a'] = in_t::A;
    ret['C'] = in_t::C;
    ret['c'] = in_t::C;
    ret['G'] = in_t::G;
    ret['g'] = in_t::G;
    ret['T'] = in_t::T;
    ret['t'] = in_t::T;
    ret['U'] = in_t::U;
    ret['u'] = in_t::U;

    // iupac characters are implicitly "UNKNOWN"
    return ret;
}()
};

//!\copydoc seqan3::dna4::complement_table
static const std::array<dna5, value_size> complement_table;

public:
//!\privatesection
//!\brief The data member.
internal_type _value;
//!\publicsection
}
;

constexpr dna5 dna5::A{ internal_type::A };
constexpr dna5 dna5::C{ internal_type::C };
constexpr dna5 dna5::G{ internal_type::G };
constexpr dna5 dna5::T{ internal_type::T };
constexpr dna5 dna5::N{ internal_type::N };
constexpr dna5 dna5::U{ dna5::T };
constexpr dna5 dna5::UNKNOWN{ dna5::N };

constexpr std::array<dna5, dna5::value_size> dna5::complement_table{
    dna5::T, // complement of dna5::A
    dna5::G, // complement of dna5::C
    dna5::C, // complement of dna5::G
    dna5::A, // complement of dna5::T
    dna5::N  // complement of dna5::N
};

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::dna5.
//!\relates dna5
using dna5_vector = std::vector<dna5>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief dna5 literal
 * \relates seqan3::dna5
 * \returns seqan3::dna5_vector
 *
 * You can use this string literal to easily assign to dna5_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna5_vector foo{"ACGTTA"};
 *     // dna5_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     dna5_vector foo{"ACGTTA"_dna5};
 *     dna5_vector bar = "ACGTTA"_dna5;
 *     auto bax = "ACGTTA"_dna5;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline dna5_vector operator""_dna5(const char * s, std::size_t n)
{
    dna5_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i) r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
