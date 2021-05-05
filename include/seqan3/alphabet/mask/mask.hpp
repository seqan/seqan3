// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
 * \ingroup mask
 * \implements seqan3::writable_semialphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
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
     * \deprecated Please use seqan3::mask::unmasked
     */
    SEQAN3_DEPRECATED_310 static const mask UNMASKED;

    /*!\brief Member for masked.
     * \details
     * \deprecated Please use seqan3::mask::masked
     */
    SEQAN3_DEPRECATED_310 static const mask MASKED;

    /*!\brief Member for unmasked.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static const mask unmasked;

    /*!\brief Member for masked.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static const mask masked;
    //!\}
};

mask constexpr mask::UNMASKED{mask{}.assign_rank(0)};
mask constexpr mask::MASKED{mask{}.assign_rank(1)};
mask constexpr mask::unmasked{mask{}.assign_rank(0)};
mask constexpr mask::masked{mask{}.assign_rank(1)};
} // namespace seqan3
