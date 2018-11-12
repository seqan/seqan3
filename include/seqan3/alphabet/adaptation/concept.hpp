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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::CharAdaptation and seqan3::UintAdaptation.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{
/*!\interface seqan3::CharAdaptation <>
 * \extends seqan3::Alphabet
 * \brief A concept that covers char type adaptations for seqan3::Alphabet.
 * \ingroup adaptation
 *
 * \details
 * This concept introduces no formal requirements beyond those of seqan3::Alphabet
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
concept CharAdaptation = Alphabet<type> &&
                                       detail::is_char_adaptation_v<std::remove_reference_t<type>>;
//!\endcond

/*!\interface seqan3::UintAdaptation <>
 * \extends seqan3::Alphabet
 * \brief A concept that covers uint type adaptations for seqan3::Alphabet.
 * \ingroup adaptation
 *
 * \details
 * This concept introduces no formal requirements beyond those of seqan3::Alphabet
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
concept UintAdaptation = Alphabet<type> &&
                                       detail::is_uint_adaptation_v<std::remove_reference_t<type>>;
//!\endcond

} // namespace seqan3
