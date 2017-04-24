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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// Author: David Heller <david.heller@fu-berlin.de>
// ==========================================================================

#pragma once

#include <optional>
#include <cassert>

#include <seqan3/alphabet/concept.hpp>

/*! The gap alphabet
 * \ingroup alphabet
 */

namespace seqan3
{

/*! The gap alphabet of -
 *
 * The alphabet always has the same value ('-').
 *
 *     gap my_gap{};
 *     gap another_gap{}.assign_char('A'); // setting this does not change anything
 *
 *     if (my_gap.to_char() == another_gap.to_char())
 *        std::cout << "Both gaps are the same!";
 */

struct gap
{
    /* types */
    //! the type of the alphabet when converted to char (e.g. via @link to_char @endlink)
    using char_type = char;
    //! the type of the alphabet when represented as a number (e.g. via @link to_rank @endlink)
    using rank_type = bool;

    /* member */
    //! internal value
    static constexpr rank_type value = 0;

    //! The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{1};

    /* public member functions */
    //! ability to cast to @link char_type @endlink **explicitly**.
    explicit constexpr operator char_type() const
    {
        return to_char();
    }

    //! return the letter as a character of @link char_type @endlink.
    constexpr char_type to_char() const
    {
        return '-';
    }

    //! return the gap's numeric value.
    constexpr rank_type to_rank() const
    {
        return value;
    }

    //! assign from a character does not do anything.
    constexpr gap & assign_char(char_type const in)
    {
        return *this;
    }

    //! assign from a numeric value does not do anything.
    constexpr gap & assign_rank(rank_type const in)
    {
        assert(value == 0);
        return *this;
    }

    //! @name comparison operators
    //!@{
    constexpr bool operator==(gap const & rhs) const
    {
        return value == rhs.value;
    }

    constexpr bool operator!=(gap const & rhs) const
    {
        return value != rhs.value;
    }

    constexpr bool operator<(gap const & rhs) const
    {
        return value < rhs.value;
    }

    constexpr bool operator>(gap const & rhs) const
    {
        return value > rhs.value;
    }

    constexpr bool operator<=(gap const & rhs) const
    {
        return value <= rhs.value;
    }

    constexpr bool operator>=(gap const & rhs) const
    {
        return value >= rhs.value;
    }
    //!@}
};

#ifndef NDEBUG
static_assert(alphabet_concept<gap>);
#endif

}
