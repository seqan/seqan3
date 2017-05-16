// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

/*!\file
 * \ingroup alphabet
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Contains seqan3::gap.
 */

#pragma once

#include <cassert>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

/*!\brief The alphabet of a gap character '-'
 *
 * The alphabet always has the same value ('-').
 *
 * ```cpp
 *     gap my_gap = gap::GAP;
 *     gap another_gap{}.assign_char('A'); // setting this does not change anything
 *
 *     std::cout << my_gap.to_char(); // outputs '-'
 *     if (my_gap.to_char() == another_gap.to_char())
 *        std::cout << "Both gaps are the same!";
 * ```
 */

struct gap
{
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = bool;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     */
    //!\{
    static const gap GAP;
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type (returns always '-').
    constexpr char_type to_char() const
    {
        return '-';
    }

    //!\brief Return the letter's numeric value or rank in the alphabet. (returns always 0)
    constexpr rank_type to_rank() const
    {
        return 0;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character (no-op, since gap has only one character).
    constexpr gap & assign_char(char_type const in)
    {
        return *this;
    }

    //!\brief Assign from a numeric value (no-op, since gap has only one character).
    constexpr gap & assign_rank(rank_type const in)
    {
        assert(in == 0);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{1};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(gap const & rhs) const
    {
        return true;
    }

    constexpr bool operator!=(gap const & rhs) const
    {
        return false;
    }

    constexpr bool operator<(gap const & rhs) const
    {
        return false;
    }

    constexpr bool operator>(gap const & rhs) const
    {
        return false;
    }

    constexpr bool operator<=(gap const & rhs) const
    {
        return true;
    }

    constexpr bool operator>=(gap const & rhs) const
    {
        return true;
    }
    //!\}
};

constexpr gap gap::GAP{};

#ifndef NDEBUG
static_assert(alphabet_concept<gap>);
#endif

}
