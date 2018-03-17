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
 * The phred42 alphabet structure represents the zero-based phred score range
 * [0..41], mapped to the ascii-ordered range ['!' .. 'J']. It therefore can
 * represent the Illumina 1.8 standard score and the original Sanger score.
*/
struct phred42
{
    /*!\name Member types
    * \{
    */
    //!\brief The 0-based integer representation of a quality score.
    using rank_type = uint8_t;

    //!\brief The integer representation of a quality score.
    using phred_type = uint8_t;

    //!\brief The '!'-based character representation of a quality score.
    using char_type = char;
    //!\}

    //!\privatesection
    /*!\name Member variables.
    * \{
    */
    //!\brief The internal 0-based rank value.
    rank_type _value;

    //!\brief The projection offset between char and rank quality score representation.
    static constexpr char_type offset_char{'!'};
    //!\}

public:
    //!\brief Value assignment with implicit compatibility to inner type.
    constexpr phred42 & operator =(phred_type const c)
    {
        assert(c < max_value_size);
        _value = (c < value_size) ? c : value_size - 1;
        return *this;
    }

    /*!\name Comparison operators.
    * \{
    */
    constexpr bool operator==(const phred42 & rhs) const
    {
        return this->_value == rhs._value;
    }

    constexpr bool operator!=(const phred42 & rhs) const
    {
        return this->_value != rhs._value;
    }

    constexpr bool operator<(const phred42 & rhs) const
    {
        return this->_value < rhs._value;
    }

    constexpr bool operator>(const phred42 & rhs) const
    {
        return this->_value > rhs._value;
    }

    constexpr bool operator<=(const phred42 & rhs) const
    {
        return this->_value <= rhs._value;
    }

    constexpr bool operator>=(const phred42 & rhs) const
    {
        return this->_value >= rhs._value;
    }
    //!\}

    //!\brief Convert quality score to its ascii representation.
    constexpr char_type to_char() const
    {
        return static_cast<char_type>(_value + offset_char);
    }

    //!\brief Set internal value given its ascii representation.
    constexpr phred42 & assign_char(char_type const c)
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Explicit compatibility to internal rank representation.
    constexpr rank_type to_rank() const
    {
        return _value;
    }

    //!\brief Set internal value given phred integer code p.
    constexpr phred42 & assign_phred(phred_type const p)
    {
        // p >= 0 is implicitly always true
        assert(p < max_value_size);
        _value = (p < value_size) ? p : value_size - 1;
        return *this;
    }

    //!\brief Set internal value given 0-based rank code p.
    constexpr phred42 & assign_rank(phred_type const p)
    {
        return assign_phred(p);
    }

    //!\brief Return the integer rank code.
    constexpr rank_type to_phred() const
    {
        return _value;
    }
    //!\}

    //!\brief The phred score range size for Illumina 1.8 standard.
    static constexpr rank_type value_size{42};

protected:
    //!\privatesection
    /*!\name Member variables.
    * \{
    */
    //!\brief Maximally tolerated phred score [42..61] that will be mapped to 41.
    static constexpr rank_type max_value_size{62};
    //!\}

    //!\brief Char to value conversion table.
    static constexpr std::array<char_type, 256> char_to_value
    {
        [] () constexpr
        {
            std::array<char_type, 256> ret{};
            for (char_type c = '!'; c <= 'J'; ++c)
                ret[c] = c - '!';
            // reduce ['K' .. '_']  to 'J'
            for (char_type c = 'K'; c <= '^'; ++c)
                ret[c] = ret['J'];
            return ret;
        }()
    };
};

} // namespace seqan3
