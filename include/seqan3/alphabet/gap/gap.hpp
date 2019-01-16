// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Contains seqan3::gap.
 */

#pragma once

#include <cassert>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The alphabet of a gap character '-'
 * \ingroup gap
 * \implements seqan3::alphabet_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
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
    constexpr gap & assign_rank([[maybe_unused]] rank_type const i) noexcept
    {
        assert(i == 0);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{1};

    //!\name Comparison operators
    //!\{
    friend constexpr bool operator==(gap const &, gap const &) noexcept
    {
        return true;
    }

    friend constexpr bool operator!=(gap const &, gap const &) noexcept
    {
        return false;
    }

    friend constexpr bool operator<(gap const &, gap const &) noexcept
    {
        return false;
    }

    friend constexpr bool operator>(gap const &, gap const &) noexcept
    {
        return false;
    }

    friend constexpr bool operator<=(gap const &, gap const &) noexcept
    {
        return true;
    }

    friend constexpr bool operator>=(gap const &, gap const &) noexcept
    {
        return true;
    }
    //!\}
};
}
