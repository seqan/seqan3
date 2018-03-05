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
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains seqan3::dna15, container aliases and string literals.
 */

#pragma once

#include <cassert>

#include <string>
#include <vector>

#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

// ------------------------------------------------------------------
// dna15
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The 15 letter DNA alphabet, containing all IUPAC smybols minus the gap.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_concept
 *
 * \details
 * Note that you can assign 'U' as a character to dna15 and it will silently
 * be converted to 'T'.
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     dna15 my_letter{dna15::A};
 *     // doesn't work:
 *     // dna15 my_letter{'A'};
 *
 *     my_letter.assign_char('C'); // <- this does!
 *
 *     my_letter.assign_char('F'); // converted to N internally
 *     if (my_letter.to_char() == 'N')
 *        std::cout << "yeah\n"; // "yeah";
 *~~~~~~~~~~~~~~~
 */

struct dna15
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
    static const dna15 A;
    static const dna15 B;
    static const dna15 C;
    static const dna15 D;
    static const dna15 G;
    static const dna15 H;
    static const dna15 K;
    static const dna15 M;
    static const dna15 N;
    static const dna15 R;
    static const dna15 S;
    static const dna15 T;
    static const dna15 U;
    static const dna15 V;
    static const dna15 W;
    static const dna15 Y;
    static const dna15 UNKNOWN;
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\copydoc seqan3::dna4::to_char
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    //!\copydoc seqan3::dna4::to_rank
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }

    //!\copydoc seqan3::dna4::complement
    constexpr dna15 complement() const noexcept
    {
        return complement_table[to_rank()];
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\copydoc seqan3::dna4::assign_char
    constexpr dna15 & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\copydoc seqan3::dna4::assign_rank
    constexpr dna15 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\copydoc seqan3::dna4::value_size
    static constexpr rank_type value_size{15};

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
    constexpr bool operator==(dna15 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(dna15 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(dna15 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(dna15 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(dna15 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(dna15 const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}

protected:
    //!\privatesection
    //!\copydoc seqan3::dna4::internal_type
    enum struct internal_type : rank_type
    {
        A,
        B,
        C,
        D,
        G,
        H,
        K,
        M,
        N,
        R,
        S,
        T,
        U = T,
        V,
        W,
        Y,
        UNKNOWN = N
    };

    //!\copydoc seqan3::dna4::value_to_char
    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'B',
        'C',
        'D',
        'G',
        'H',
        'K',
        'M',
        'N',
        'R',
        'S',
        'T',
        'V',
        'W',
        'Y'
    };

    //!\copydoc seqan3::dna4::char_to_value
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

            // iupac characters
            ret['R'] = in_t::R; ret['r'] = in_t::R;
            ret['Y'] = in_t::Y; ret['y'] = in_t::Y;
            ret['S'] = in_t::S; ret['s'] = in_t::S;
            ret['W'] = in_t::W; ret['w'] = in_t::W;
            ret['K'] = in_t::K; ret['k'] = in_t::K;
            ret['M'] = in_t::M; ret['m'] = in_t::M;
            ret['B'] = in_t::B; ret['b'] = in_t::B;
            ret['D'] = in_t::D; ret['d'] = in_t::D;
            ret['H'] = in_t::H; ret['h'] = in_t::H;
            ret['V'] = in_t::V; ret['v'] = in_t::V;
            ret['N'] = in_t::N; ret['n'] = in_t::N;
            return ret;
        }()
    };

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<dna15, value_size> complement_table;

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr dna15 dna15::A{internal_type::A};
constexpr dna15 dna15::B{internal_type::B};
constexpr dna15 dna15::C{internal_type::C};
constexpr dna15 dna15::D{internal_type::D};
constexpr dna15 dna15::G{internal_type::G};
constexpr dna15 dna15::H{internal_type::H};
constexpr dna15 dna15::K{internal_type::K};
constexpr dna15 dna15::M{internal_type::M};
constexpr dna15 dna15::N{internal_type::N};
constexpr dna15 dna15::R{internal_type::R};
constexpr dna15 dna15::S{internal_type::S};
constexpr dna15 dna15::T{internal_type::T};
constexpr dna15 dna15::U{dna15::T};
constexpr dna15 dna15::V{internal_type::V};
constexpr dna15 dna15::W{internal_type::W};
constexpr dna15 dna15::Y{internal_type::Y};
constexpr dna15 dna15::UNKNOWN{dna15::N};

constexpr std::array<dna15, dna15::value_size> dna15::complement_table
{
    dna15::T,    // complement of dna15::A
    dna15::V,    // complement of dna15::B
    dna15::G,    // complement of dna15::C
    dna15::H,    // complement of dna15::D
    dna15::C,    // complement of dna15::G
    dna15::D,    // complement of dna15::H
    dna15::M,    // complement of dna15::K
    dna15::K,    // complement of dna15::M
    dna15::N,    // complement of dna15::N
    dna15::Y,    // complement of dna15::R
    dna15::S,    // complement of dna15::S
    dna15::A,    // complement of dna15::T
    dna15::B,    // complement of dna15::V
    dna15::W,    // complement of dna15::W
    dna15::R     // complement of dna15::Y
};

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::dna15.
//!\relates dna15
using dna15_vector = std::vector<dna15>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief dna15 literal
 * \relates seqan3::dna15
 * \returns seqan3::dna15_vector
 *
 * You can use this string literal to easily assign to dna15_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     // these don't work:
 *     // dna15_vector foo{"ACGTTA"};
 *     // dna15_vector bar = "ACGTTA";
 *
 *     // but these do:
 *     using namespace seqan3::literal;
 *     dna15_vector foo{"ACGTTA"_dna15};
 *     dna15_vector bar = "ACGTTA"_dna15;
 *     auto bax = "ACGTTA"_dna15;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline dna15_vector operator""_dna15(const char * s, std::size_t n)
{
    dna15_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
