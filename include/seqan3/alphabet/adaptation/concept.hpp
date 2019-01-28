// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::char_adaptation_concept and seqan3::uint_adaptation_concept.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{
/*!\interface seqan3::char_adaptation_concept <>
 * \extends seqan3::alphabet_concept
 * \brief A concept that covers char type adaptations for seqan3::alphabet_concept.
 * \ingroup adaptation
 *
 * \details
 * This concept introduces no formal requirements beyond those of seqan3::alphabet_concept
 * and type being one of the following types:
 *
 *   * `char`
 *   * `char16_t`
 *   * `char32_t`
 *
 * \attention
 *   * Note that `signed char` and `unsigned char` are absent from the list, because of
 * their type ambiguity with `int8_t` and `uint8_t`.
 *   * Note that `wchar_t` is absent from the list for its notorious brokenness (different sizes and signedness
 * between platforms); use `char16_t` or `char32_t` instead.
 *
 * \attention
 * Please be aware that this file needs be included **after** `alphabet/adaptation/char.hpp`.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT char_adaptation_concept = alphabet_concept<type> &&
                                  detail::is_char_adaptation_v<type>;
//!\endcond

/*!\interface seqan3::uint_adaptation_concept <>
 * \extends seqan3::alphabet_concept
 * \brief A concept that covers uint type adaptations for seqan3::alphabet_concept.
 * \ingroup adaptation
 *
 * \details
 * This concept introduces no formal requirements beyond those of seqan3::alphabet_concept
 * and type being one of the following types:
 *
 *   * `uint8_t`
 *   * `uint16_t`
 *   * `uint32_t`
 *
 * \attention
 * Note that `uint64_t` is absent from the list, because there is no corresponding
 * character type.
 *
 * \attention
 * Please be aware that this file needs be included **after** `alphabet/adaptation/uint.hpp`.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT uint_adaptation_concept = alphabet_concept<type> &&
                                  detail::is_uint_adaptation_v<type>;
//!\endcond

} // namespace seqan3
