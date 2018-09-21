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
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Contains seqan3::aa20, container aliases and string literals.
 */

#pragma once

#include <cassert>
#include <vector>

#include <seqan3/core/platform.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>

namespace seqan3
{
/*!\brief The canonical amino acid alphabet.
 * \ingroup aminoacid
 * \implements seqan3::aminoacid_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 *
 * \details
 * The alphabet consists of letters A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
 *
 * The alphabet may be brace initialized from the static letter members (see above). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 * *Note:* Letters which belong in the extended alphabet will be automatically converted based on the frequency
 * of their options.\n Terminator characters are converted to W, because the most commonly occurring
 * stop codon in higher eukaryotes is UGA<sup>2</sup>.
 * Anything unknown is converted to S, because it occurs most frequently across 53 vertebrates<sup>1</sup>.
 *
 * |Input Letter  |Converts to  |
 * |--------------|-------------|
 * |B             |D<sup>1</sup>|
 * |J             |L<sup>1</sup>|
 * |O             |L<sup>1</sup>|
 * |U             |C<sup>1</sup>|
 * |Z             |E<sup>1</sup>|
 * |X (Unknown)   |S<sup>1</sup>|
 * |* (Terminator)|W<sup>2</sup>|
 * <sup><b>1</b></sup>King, J. L., & Jukes, T. H. (1969). Non-Darwinian Evolution.
 * Science, 164(3881), 788-798. doi:10.1126/science.164.3881.788\n
 * <sup><b>2</b></sup>Trotta, E. (2016). Selective forces and mutational biases drive stop codon usage
 * in the human genome: a comparison with sense codon usage.
 * BMC Genomics, 17, 366. http://doi.org/10.1186/s12864-016-2692-4
 *
 * \snippet test/snippet/alphabet/aminoacid/aa20.cpp construction
 */

struct aa20
{
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;

    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = uint8_t;


    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.
     */
    //!\{
    static const aa20 A;
    static const aa20 B;
    static const aa20 C;
    static const aa20 D;
    static const aa20 E;
    static const aa20 F;
    static const aa20 G;
    static const aa20 H;
    static const aa20 J;
    static const aa20 I;
    static const aa20 K;
    static const aa20 L;
    static const aa20 M;
    static const aa20 N;
    static const aa20 O;
    static const aa20 P;
    static const aa20 Q;
    static const aa20 R;
    static const aa20 S;
    static const aa20 T;
    static const aa20 U;
    static const aa20 V;
    static const aa20 W;
    static const aa20 X;
    static const aa20 Y;
    static const aa20 Z;
    static const aa20 TERMINATOR;
    static const aa20 UNKNOWN;
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
    constexpr aa20 & assign_char(char_type const c) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(c)];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr aa20 & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{20};

    /*!\name Conversion operators
     * \{
     */
    //!\brief Explicit conversion to any other amino acid alphabet (via char representation).
    //!\tparam other_aa_type The type to convert to; must satisfy seqan3::aminoacid_concept.
    template <typename other_aa_type>
    //!\cond
        requires aminoacid_concept<other_aa_type>
    //!\endcond
    explicit constexpr operator other_aa_type() const noexcept
    {
        return detail::convert_through_char_representation<other_aa_type, std::decay_t<decltype(*this)>>[to_rank()];
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(aa20 const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(aa20 const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(aa20 const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(aa20 const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(aa20 const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(aa20 const & rhs) const noexcept
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
     * the static "letter values" members necessary.
     */
    enum struct internal_type : rank_type
    {
        A,
        C,
        D,
        E,
        F,
        G,
        H,
        I,
        K,
        L,
        M,
        N,
        P,
        Q,
        R,
        S,
        T,
        V,
        W,
        Y,
        B = D,
        J = L,
        O = L,
        U = C,
        X = S,
        Z = E,
        TERMINATOR = W,
        UNKNOWN = S
    };

    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        'A',
        'C',
        'D',
        'E',
        'F',
        'G',
        'H',
        'I',
        'K',
        'L',
        'M',
        'N',
        'P',
        'Q',
        'R',
        'S',
        'T',
        'V',
        'W',
        'Y',
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

            ret['A'] = in_t::A; ret['a'] = in_t::A;
            ret['B'] = in_t::D; ret['b'] = in_t::D; // Convert b (either D/N) to D, since D occurs more frequently.
            ret['C'] = in_t::C; ret['c'] = in_t::C;
            ret['D'] = in_t::D; ret['d'] = in_t::D;
            ret['E'] = in_t::E; ret['e'] = in_t::E;
            ret['F'] = in_t::F; ret['f'] = in_t::F;
            ret['G'] = in_t::G; ret['g'] = in_t::G;
            ret['H'] = in_t::H; ret['h'] = in_t::H;
            ret['I'] = in_t::I; ret['i'] = in_t::I;
            ret['J'] = in_t::L; ret['j'] = in_t::L; // Convert j (either I/L) to L, since L occurs more frequently.
            ret['K'] = in_t::K; ret['k'] = in_t::K;
            ret['L'] = in_t::L; ret['l'] = in_t::L;
            ret['M'] = in_t::M; ret['m'] = in_t::M;
            ret['N'] = in_t::N; ret['n'] = in_t::N;
            ret['O'] = in_t::L; ret['o'] = in_t::L; // Convert Pyrrolysine to lysine.
            ret['P'] = in_t::P; ret['p'] = in_t::P;
            ret['Q'] = in_t::Q; ret['q'] = in_t::Q;
            ret['R'] = in_t::R; ret['r'] = in_t::R;
            ret['S'] = in_t::S; ret['s'] = in_t::S;
            ret['T'] = in_t::T; ret['t'] = in_t::T;
            ret['U'] = in_t::C; ret['u'] = in_t::C; // Convert Selenocysteine to cysteine.
            ret['V'] = in_t::V; ret['v'] = in_t::V;
            ret['W'] = in_t::W; ret['w'] = in_t::W;
            ret['X'] = in_t::S; ret['x'] = in_t::S; // Convert unknown amino acids to serine.
            ret['Y'] = in_t::Y; ret['y'] = in_t::Y;
            ret['Z'] = in_t::E; ret['z'] = in_t::E; // Convert z (either E/Q) to E, since E occurs more frequently.
            ret['*'] = in_t::W; // The most common stop codon is UGA. This is most similar to a Tryptophan.
            return ret;
        }()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr aa20 aa20::A{internal_type::A};
constexpr aa20 aa20::C{internal_type::C};
constexpr aa20 aa20::D{internal_type::D};
constexpr aa20 aa20::E{internal_type::E};
constexpr aa20 aa20::F{internal_type::F};
constexpr aa20 aa20::G{internal_type::G};
constexpr aa20 aa20::H{internal_type::H};
constexpr aa20 aa20::I{internal_type::I};
constexpr aa20 aa20::K{internal_type::K};
constexpr aa20 aa20::L{internal_type::L};
constexpr aa20 aa20::M{internal_type::M};
constexpr aa20 aa20::N{internal_type::N};
constexpr aa20 aa20::P{internal_type::P};
constexpr aa20 aa20::Q{internal_type::Q};
constexpr aa20 aa20::R{internal_type::R};
constexpr aa20 aa20::S{internal_type::S};
constexpr aa20 aa20::T{internal_type::T};
constexpr aa20 aa20::V{internal_type::V};
constexpr aa20 aa20::W{internal_type::W};
constexpr aa20 aa20::Y{internal_type::Y};
constexpr aa20 aa20::B{aa20::D};
constexpr aa20 aa20::J{aa20::L};
constexpr aa20 aa20::O{aa20::L};
constexpr aa20 aa20::U{aa20::C};
constexpr aa20 aa20::X{aa20::S};
constexpr aa20 aa20::Z{aa20::E};
constexpr aa20 aa20::TERMINATOR{aa20::W};
constexpr aa20 aa20::UNKNOWN{aa20::S};

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{
//!\brief Alias for an std::vector of seqan3::aa20.
//!\relates aa20
using aa20_vector = std::vector<aa20>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief aa20 literal
 * \relates seqan3::aa20
 * \returns seqan3::aa20_vector
 *
 * You can use this string literal to easily assign to aa20_vector:
 *
 *\snippet test/snippet/alphabet/aminoacid/aa20.cpp literal
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline aa20_vector operator""_aa20(const char * s, std::size_t n)
{
    aa20_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
