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

#include <seqan3/alphabet/detail/alphabet_base.hpp>

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

class gap : public alphabet_base<gap, 1, char>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<gap, 1, char>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

    //!\brief The character that will be printed.
    static constexpr char char_value = '-';

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr gap() noexcept : base_t{} {}
    constexpr gap(gap const &) = default;
    constexpr gap(gap &&) = default;
    constexpr gap & operator=(gap const &) = default;
    constexpr gap & operator=(gap &&) = default;
    ~gap() = default;

    using base_t::base_t;
    //!\}
};

}
