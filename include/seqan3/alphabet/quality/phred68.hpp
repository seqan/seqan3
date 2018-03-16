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
// Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================

#pragma once

#include <cassert>

namespace seqan3
{

/*!
 * The phred68 alphabet represents the phred score range [-5..62] mapped
 * to the ascii-ordered range [';' .. '~']. It represents the Solexa standard
 * score.
*/
struct phred68
{
    //!\publicsection
    /*!\name Member types
    * \{
    */
    //!\brief The 0-based integer representation of a quality score.
    //!\hideinitializer
    using rank_type = uint8_t;

    //!\brief The integer representation of a quality score.
    //!\hideinitializer
    using phred_type = int8_t;

    //!\brief The ';'-based character representation of a quality score.
    //!\hideinitializer
    using char_type = char;
    //!\}

    //!\privatesection
    /*!\name Member variables.
    * \{
    */
    //!\brief The internal 0-based rank value.
    rank_type _value;

    //!\brief The projection offset between char and rank quality score representation.
    //!\hideinitializer
    static constexpr char_type offset_char{';'};

    //! projection offsets of a phred quality score
    static constexpr phred_type offset_phred{-5};
    //!\}

    //!\brief Value assignment with implicit compatibility to inner type.
    constexpr phred68 & operator =(phred_type const c)
    {
        assert((c >= offset_phred) & (c < value_size + offset_phred));
        _value = c - offset_phred;  // assign its rank value
        return *this;
    }

    /*!\name Comparison operators.
    * \{
    */
    constexpr bool operator==(const phred68 & rhs) const
    {
        return this->_value == rhs._value;
    }

    constexpr bool operator!=(const phred68 & rhs) const
    {
        return this->_value != rhs._value;
    }

    constexpr bool operator<(const phred68 & rhs) const
    {
        return this->_value < rhs._value;
    }

    constexpr bool operator>(const phred68 & rhs) const
    {
        return this->_value > rhs._value;
    }

    constexpr bool operator<=(const phred68 & rhs) const
    {
        return this->_value <= rhs._value;
    }

    constexpr bool operator>=(const phred68 & rhs) const
    {
        return this->_value >= rhs._value;
    }
    //!\}

    /*!\name Conversion and explicit conversion functions.
    * \{
    */
    //!\brief Convert quality score to its ascii representation.
    constexpr char_type to_char() const
    {
        return static_cast<char_type>(_value + offset_char);
    }

    //!\brief Set internal value given its ascii representation.
    constexpr phred68 & assign_char(char_type const c)
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Explicit compatibility to internal rank representation.
    constexpr rank_type to_rank() const
    {
        return _value;
    }

    //!\brief Set internal value given phred_offset-based integer code p.
    constexpr phred68 & assign_phred(phred_type const p)
    {
        // p >= 0 is implicitly always true
        assert((p >= offset_phred) & (p < value_size + offset_phred));
        _value = p - offset_phred;
        return *this;
    }

    //!\brief Set internal value given 0-based rank code p.
    constexpr phred68 & assign_rank(rank_type const p)
    {
        // p implicitly always true due to unsignedness of p
        assert(p < value_size);
        _value = p;
        return *this;
    }

    //!\brief Return the integer phred score [offset_phred .. offset_phred+68[.
    constexpr rank_type to_phred() const
    {
        return _value + offset_phred;
    }
    //!\}

    //!\brief The phred score range size for the Solexa standard.
    static constexpr rank_type value_size{68};

protected:
    //!\brief Char to value conversion table.
    static constexpr std::array<char_type, 256> char_to_value
    {
        [] () constexpr
        {
            std::array<char_type, 256> ret{};
            for (char_type c = ';'; c <= '~'; ++c)
                ret[c] = c - offset_char;
            return ret;
        }()
    };
};

} // namespace seqan3
