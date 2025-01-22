// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Create a mask composite which can be applied with another alphabet.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 */

#pragma once

#include <cassert>

#include <seqan3/alphabet/alphabet_base.hpp>

namespace seqan3
{
/*!\brief Implementation of a masked alphabet to be used for tuple composites.
 * \ingroup alphabet_mask
 * \implements seqan3::writable_semialphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 *
 * \details
 *
 * This alphabet is not usually used directly, but instead via seqan3::masked.
 *
 * \include test/snippet/alphabet/mask/mask.cpp
 *
 * \see \link mask Mask submodule \endlink, it contains an explanation of hard-masking (unknown character) and
 *      soft-masking (lower/upper case letters).
 *
 * \stableapi{Since version 3.1.}
 */
class mask : public alphabet_base<mask, 2, void>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<mask, 2, void>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr mask() = default;                         //!< Defaulted.
    constexpr mask(mask const &) = default;             //!< Defaulted.
    constexpr mask(mask &&) = default;                  //!< Defaulted.
    constexpr mask & operator=(mask const &) = default; //!< Defaulted.
    constexpr mask & operator=(mask &&) = default;      //!< Defaulted.
    ~mask() = default;                                  //!< Defaulted.

    //!\}

    /*!\name Boolean values
     * \brief Static member "booleans" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface.
     */
    //!\{
    /*!\brief Member for unmasked.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static mask const unmasked;

    /*!\brief Member for masked.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static mask const masked;
    //!\}
};

constexpr mask mask::unmasked{mask{}.assign_rank(0)};
constexpr mask mask::masked{mask{}.assign_rank(1)};
} // namespace seqan3
