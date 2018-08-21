// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

/*!\file
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
 * \ingroup gap
 * \implements seqan3::alphabet_concept
 *
 * The alphabet always has the same value ('-').
 *
 * \snippet test/snippet/alphabet/gap/gap.cpp general
 */

struct gap
{
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = bool;

    //!\cond
    bool _bug_workaround; // See GCC Bug-Report: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87113
    //!\endcond

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
    constexpr char_type to_char() const noexcept
    {
        return '-';
    }

    //!\brief Return the letter's numeric value or rank in the alphabet. (returns always 0)
    constexpr rank_type to_rank() const noexcept
    {
        return 0;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character (no-op, since gap has only one character).
    //!\param c not used, since gap has only one character
    constexpr gap & assign_char([[maybe_unused]] char_type const c) noexcept
    {
        return *this;
    }

    //!\brief Assign from a numeric value (no-op, since gap has only one character).
    //!\param i not used, since gap has only one character
    constexpr gap & assign_rank([[maybe_unused]] rank_type const i) /*noexcept*/
    {
        // TODO(marehr): mark function noexcept if assert is replaced
        // https://github.com/seqan/seqan3/issues/85
        assert(i == 0);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{1};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(gap const &) const noexcept
    {
        return true;
    }

    constexpr bool operator!=(gap const &) const noexcept
    {
        return false;
    }

    constexpr bool operator<(gap const &) const noexcept
    {
        return false;
    }

    constexpr bool operator>(gap const &) const noexcept
    {
        return false;
    }

    constexpr bool operator<=(gap const &) const noexcept
    {
        return true;
    }

    constexpr bool operator>=(gap const &) const noexcept
    {
        return true;
    }
    //!\}
};

constexpr gap gap::GAP{};

}
