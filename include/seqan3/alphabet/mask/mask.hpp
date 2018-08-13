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
 * \brief Create a mask composition which can be applied with another alphabet.
 */

#pragma once

#include <cassert>
#include <seqan3/alphabet/concept_pre.hpp>

namespace seqan3
{
/*!\brief Implementation of a masked alphabet to be used for cartesian compositions.
 * \ingroup mask
 * \implements seqan3::semi_alphabet_concept
 *
 * \details
 * This alphabet is not usually used directly, but instead via seqan3::masked.
 * For more information see the \link mask Mask submodule \endlink.
 */
struct mask
{
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = bool;

    /*!\name Boolean values
     * \brief Static member "booleans" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface . *Don't worry about the `internal_type`.*
     */
    //!\{
    static const mask UNMASKED;
    static const mask MASKED;
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr uint8_t value_size{2};

    /*!\name Read function
     * \{
     */
    //!\brief Return the letter's numeric value or rank in the alphabet.
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a numeric value or true/false.
    constexpr mask & assign_rank(rank_type const c) noexcept
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(mask const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(mask const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(mask const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(mask const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(mask const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(mask const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}
protected:
    //!\privatesection
    /*!\brief The internal type is a strictly typed enum.
     *
     * This is done to prevent aggregate initialization from numbers and/or chars.
     * It has the drawback that it also introduces a scope which in turn makes
     * the static "letter values" members necessary.
     */
    enum struct internal_type : rank_type
    {
        UNMASKED,
        MASKED
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

constexpr mask mask::UNMASKED{internal_type::UNMASKED};
constexpr mask mask::MASKED{internal_type::MASKED};
} // namespace seqan3
