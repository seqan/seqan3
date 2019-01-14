// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Create a mask composition which can be applied with another alphabet.
 */

#pragma once

#include <cassert>
#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/detail/alphabet_base.hpp>

namespace seqan3
{
/*!\brief Implementation of a masked alphabet to be used for cartesian compositions.
 * \ingroup mask
 * \implements seqan3::semi_alphabet_concept
 * \implements seqan3::detail::semi_constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \details
 * This alphabet is not usually used directly, but instead via seqan3::masked.
 * For more information see the \link mask Mask submodule \endlink.
 *
 * \snippet test/snippet/alphabet/mask/mask.cpp general
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
    constexpr mask() : base_t{} {}
    constexpr mask(mask const &) = default;
    constexpr mask(mask &&) = default;
    constexpr mask & operator=(mask const &) = default;
    constexpr mask & operator=(mask &&) = default;
    ~mask() = default;
    //!\}

    /*!\name Boolean values
     * \brief Static member "booleans" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface.
     */
    //!\{
    static const mask UNMASKED;
    static const mask MASKED;
    //!\}
};

mask constexpr mask::UNMASKED{mask{}.assign_rank(0)};
mask constexpr mask::MASKED  {mask{}.assign_rank(1)};
} // namespace seqan3
