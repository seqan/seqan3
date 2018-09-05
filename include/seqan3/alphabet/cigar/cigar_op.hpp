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
 * \brief Introduces the cigar_op alphabet.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <type_traits>

// ------------------------------------------------------------------
// cigar_op
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The (extended) cigar operation alphabet of M,D,I,H,N,P,S,X,EQ.
 * \ingroup cigar
 *
 * \details
 *
 * The CIGAR string can be either basic or extended. The only difference in the
 * extended cigar alphabet is that aligned bases are classified as an actual
 * match (letter: seqan3::cigar_op::EQ, char representation '=') or mismatch
 * (letter: seqan3::cigar_op::X, char 'X'). In contrast, the basic cigar alphabet
 * only indicated the aligned status with an seqan3::cigar_op::M, without further
 * information if the bases are actually equal or not.
 *
 * The main purpose of the seqan3::cigar_op alphabet is to be used in the
 * seqan3::cigar composition, where a cigar operation is paired with a count
 * value.
 *
 * Example usage:
 * \snippet test/snippet/alphabet/cigar/cigar_op.cpp general
 *
 * \note Usually you do not want to manipulate cigar elements and vectors on
 *       your own but convert an alignment to a cigar and back. See
 *       seqan3::get_cigar_vector for how to convert two aligned sequences into
 *       an cigar_vector.
 */
struct cigar_op
{
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = uint8_t;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used
     *        in aggregate initialization.
     * \details Similar to an Enum interface. *Don't worry about the `internal_type`.*
     */
    //!\{
    //!\brief Match/Aligned Bases; can be either an alignment match or mismatch of two aligned bases.
    static const cigar_op M;
    //!\brief Deletion; the nucleotide is present in the reference but not in the read
    static const cigar_op D;
    //!\brief Insertion; the nucleotide is present in the read but not in the reference.
    static const cigar_op I;
    //!\brief Soft Clipping; the clipped nucleotides are present in the read sequence but not part of the alignment.
    static const cigar_op S;
    //!\brief Hard Clipping; the clipped nucleotides are not present in the read.
    static const cigar_op H;
    //!\brief Skipped region; a region of nucleotides is not present in the read
    static const cigar_op N;
    //!\brief Padding; padded area in the read and not in the reference
    static const cigar_op P;
    //!\brief Alignment Mismatch; the aligned characters are not equal.
    static const cigar_op X;
    //!\brief Alignment Match; the aligned characters are equal.
    static const cigar_op EQ;
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the letter as a character of char_type.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::to_char() requirement via the seqan3::to_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    /*!\brief Return the letter's numeric value or rank in the alphabet.
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::to_rank() requirement via the seqan3::to_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a character.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::assign_char() requirement via the seqan3::assign_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr cigar_op & assign_char(char_type const c) noexcept
    {
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(c)];
        return *this;
    }

    /*!\brief Assign from a numeric value.
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::assign_rank() requirement via the seqan3::assign_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr cigar_op & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{9};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(cigar_op const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(cigar_op const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(cigar_op const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(cigar_op const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(cigar_op const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(cigar_op const & rhs) const noexcept
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
        M,
        D,
        I,
        H,
        N,
        P,
        S,
        X,
        EQ,
    };

    //!\brief Value to char conversion table.
    static constexpr char_type value_to_char[value_size]
    {
        'M',
        'D',
        'I',
        'H',
        'N',
        'P',
        'S',
        'X',
        '='
    };

    //!\brief Char to value conversion table.
    static constexpr std::array<internal_type, 256> char_to_value
    {
        [] () constexpr
        {
            using in_t = internal_type;
            std::array<in_t, 256> ret{};

            // initialize with Match M (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = in_t::M;

            // canonical
            ret['M'] = in_t::M; ret['m'] = in_t::M;
            ret['D'] = in_t::D; ret['d'] = in_t::D;
            ret['I'] = in_t::I; ret['i'] = in_t::I;
            ret['H'] = in_t::H; ret['h'] = in_t::H;
            ret['N'] = in_t::N; ret['n'] = in_t::N;
            ret['P'] = in_t::P; ret['p'] = in_t::P;
            ret['S'] = in_t::S; ret['s'] = in_t::S;
            ret['X'] = in_t::X; ret['x'] = in_t::X;
            ret['='] = in_t::EQ;

            return ret;
        }()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr cigar_op cigar_op::M{internal_type::M};
constexpr cigar_op cigar_op::D{internal_type::D};
constexpr cigar_op cigar_op::I{internal_type::I};
constexpr cigar_op cigar_op::H{internal_type::H};
constexpr cigar_op cigar_op::N{internal_type::N};
constexpr cigar_op cigar_op::P{internal_type::P};
constexpr cigar_op cigar_op::S{internal_type::S};
constexpr cigar_op cigar_op::X{internal_type::X};
constexpr cigar_op cigar_op::EQ{internal_type::EQ};

} // namespace seqan3
