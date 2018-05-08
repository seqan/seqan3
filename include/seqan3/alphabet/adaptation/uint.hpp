// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides alphabet adaptations for standard uint types.
 * \details
 * This file provides function and metafunction overloads so that the following types
 * fulfil the seqan3::alphabet_concept:
 *   * `uint8_t`
 *   * `uint16_t`
 *   * `uint32_t`
 *
 * You will likely not use these interfaces directly, they are, however, very helpful
 * for conversions between other alphabets and between other alphabets and characters.
 *
 * \attention
 * Note that `uint64_t` is absent from the list, because there is no corresponding
 * character type.
 */

#pragma once

#include <limits>

#include <meta/meta.hpp>

#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/alphabet/concept_pre.hpp>

namespace seqan3::detail
{
//!\addtogroup adaptation
//!\{

//!\brief The list of types that are defined as uint adaptations.
using uint_adaptations            = meta::list<uint8_t,
                                               uint16_t,
                                               uint32_t>;
//!\brief The corresponding list of rank types
using uint_adaptations_char_types = meta::list<char,
                                               char16_t,
                                               char32_t>;

//!\brief Metafunction overload for types that are in seqan3::detail::uint_adaptations.
template <typename type_in_list>
    requires meta::in<uint_adaptations, type_in_list>::value
struct is_uint_adaptation<type_in_list> :
    std::true_type
{};

//!\}
} // namespace seqan3::detail

namespace seqan3
{

//!\addtogroup adaptation
//!\{

// ------------------------------------------------------------------
// type metafunctions
// ------------------------------------------------------------------

/*!\brief Specialisation of seqan3::underlying_char for uint types.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \relates seqan3::uint_adaptation_concept
 * \sa seqan3::underlying_char_t
 */
template <typename uint_type>
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
struct underlying_char<uint_type>
{
    //!\brief The character type of the same size as `uint_type`.
    using type = meta::at<detail::uint_adaptations_char_types,
                          meta::find_index<detail::uint_adaptations, uint_type>>;
};

/*!\brief Specialisation of seqan3::underlying_rank for uint types.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \relates seqan3::uint_adaptation_concept
 * \sa seqan3::underlying_rank_t
 */
template <typename uint_type>
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
struct underlying_rank<uint_type>
{
    //!\brief The same as `uint_type`.
    using type = uint_type;
};

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

/*!\brief Specialisation of seqan3::alphabet_size that delegates for uint types.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \relates seqan3::uint_adaptation_concept
 * \sa seqan3::alphabet_size_v
 */
template <typename uint_type>
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
struct alphabet_size<uint_type>
{
    //!\brief Smallest unsigned integral type that can hold value;
    using type = detail::min_viable_uint_t<static_cast<uint64_t>(std::numeric_limits<uint_type>::max()) + 1 -
                                           std::numeric_limits<uint_type>::lowest()>;
    //!\brief The alphabet's size.
    static constexpr type value =
        static_cast<type>(std::numeric_limits<uint_type>::max()) + 1 - std::numeric_limits<uint_type>::lowest();
};

// ------------------------------------------------------------------
// free functions
// ------------------------------------------------------------------

/*!\name Free function wrappers for the uint alphabet adaptation
 * \brief For `uint8_t`, `uint16_t` and `uint32_t` do conversion to/from char types.
 * \ingroup adaptation
 * \relates seqan3::uint_adaptation_concept
 * \{
 */
/*!\brief Converting uint to char casts to a character type of same size.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `uint`).
 */
template <typename uint_type>
constexpr underlying_char_t<uint_type> to_char(uint_type const intgr)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr;
}

/*!\brief Converting uint to rank is a no-op (it will just return the value you pass in).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename uint_type>
constexpr underlying_rank_t<uint_type> to_rank(uint_type const intgr)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr;
}

/*!\brief Assign from a character type via implicit or explicit cast.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to assign to.
 * \param chr The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename uint_type>
constexpr uint_type & assign_char(uint_type & intgr, underlying_char_t<uint_type> const chr)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr = chr;
}

/*!\brief Assign from a character type via implicit or explicit cast.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr An alphabet letter temporary.
 * \param chr The `char` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from uint:
 * ```cpp
 * auto l = assign_char(dna5{}, 'G');  // l is of type dna5
 * ```
 */
template <typename uint_type>
constexpr uint_type && assign_char(uint_type && intgr, underlying_char_t<uint_type> const chr)
    requires detail::is_uint_adaptation_v<std::remove_reference_t<uint_type>>
{
    return std::move(intgr = chr);
}

/*!\brief Assign a rank to to the uint (same as calling `=`).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to assign to.
 * \param intgr2 The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename uint_type>
constexpr uint_type & assign_rank(uint_type & intgr, underlying_rank_t<uint_type> const intgr2)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr = intgr2;
}

/*!\brief Assign a rank to to the uint (same as calling `=`).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr An alphabet letter temporary.
 * \param intgr2 The `rank` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from rank:
 * ```cpp
 * auto l = assign_rank(dna5{}, 1);  // l is of type dna5 and == dna5::C
 * ```
 */
template <typename uint_type>
constexpr uint_type && assign_rank(uint_type && intgr, underlying_rank_t<uint_type> const intgr2)
    requires detail::is_uint_adaptation_v<std::remove_reference_t<uint_type>>
{
    return std::move(intgr = intgr2);
}
//!\}
//!\}

} // namespace seqan3

// concept tests in alphabet/adaptation.hpp because of include order constraints
