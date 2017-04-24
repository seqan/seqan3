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
 * Implementation of the Illumina 1.8 standard fulfilling the quality concept.
 * The permitted phred score range is [0 .. 41], mapped to ascii-ordered range ['!' .. 'J'].
 * For this standard internal and rank phred representation are both zero-based.
 */
struct illumina18
{
    //! the 3 representation types of a quality score
    using phred_type = int8_t;
    using rank_type = uint8_t;
    using char_type = char;

    //! internal rank value representation
    rank_type value;

    //! projection offsets of char and rank quality score
    static constexpr char_type offset_char{'!'};
    static constexpr phred_type offset_phred{0};

    //! implicit compatibility to inner_type
    constexpr illumina18 & operator =(rank_type const c)
    {
        value = c;
        return *this;
    }

    //! comparison operators
    constexpr bool operator==(const illumina18 & rhs) const
    {
        return this->value == rhs.value;
    }

    constexpr bool operator!=(const illumina18 & rhs) const
    {
        return this->value != rhs.value;
    }

    constexpr bool operator<(const illumina18 & rhs) const
    {
        return this->value < rhs.value;
    }

    constexpr bool operator>(const illumina18 & rhs) const
    {
        return this->value > rhs.value;
    }

    constexpr bool operator<=(const illumina18 & rhs) const
    {
        return this->value <= rhs.value;
    }

    constexpr bool operator>=(const illumina18 & rhs) const
    {
        return this->value >= rhs.value;
    }

    //! explicit compatibility to char code of a quality score
    explicit constexpr operator char() const
    {
        return to_char();
    }

    //! convert quality score to its 1-letter code
    constexpr char_type to_char() const
    {
        return value + offset_char;
    }

    //! set internal value given 1-letter code
    constexpr illumina18 & assign_char(char_type const c)
    {
        assert(c >= '!' && c <= 'J');
        value = c - '!';
        return *this;
    }

    //! explicit compatibility to internal rank representation
    constexpr rank_type to_rank() const
    {
        return value;
    }

    //! set internal value given zero-based integer c
    constexpr illumina18 & assign_rank(rank_type const c)
    {
        assert(c >= 0 && c < value_size);
        value = c;
        return *this;
    }

    //! set internal value given Illumina 1.8 integer code p
    constexpr illumina18 & assign_phred(phred_type const p)
    {
        assert(p >= offset_phred && p < offset_phred + value_size);
        value = p - offset_phred;
        return *this;
    }

    //! get Illumina 1.8 integer code
    constexpr phred_type to_phred() const
    {
        return value + offset_phred;
    }

    //! phred score range for Illumina 1.8 standard
    static constexpr rank_type value_size{42};
};

} // namespace seqan3

#ifndef NDEBUG
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
static_assert(seqan3::alphabet_concept<seqan3::illumina18>);
static_assert(seqan3::quality_concept<seqan3::illumina18>);
static_assert(seqan3::detail::internal_quality_concept<seqan3::illumina18>);
#endif
