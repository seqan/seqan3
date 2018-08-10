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
 * \brief Provides alphabet adaptations for standard char types.
 * \details
 * This file provides function and metafunction overloads so that the following types
 * fulfil the seqan3::alphabet_concept:
 *   * `char`
 *   * `char16_t`
 *   * `char32_t`
 *
 * You will likely not use these interfaces directly, they are, however, very helpful
 * for conversions between other alphabets and between other alphabets and characters.
 *
 * \attention
 *   * Note that `signed char` and `unsigned char` are absent from the list, because of
 * their type ambiguity with `int8_t` and `uint8_t`.
 */

#pragma once

#include <limits>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/core/detail/int_types.hpp>

namespace seqan3::detail
{
//!\addtogroup adaptation
//!\{

//!\brief The list of types that are defined as char adaptations.
using char_adaptations            = meta::list<char,
                                               char16_t,
                                               char32_t,
                                               wchar_t>;
// TODO Windows adaption requires different handling of wchar_t, in Windows it is 16 bits and holds UTF-16 code units
//!\brief The corresponding list of rank types
using char_adaptations_rank_types = meta::list<std::uint_least8_t,
                                               std::uint_least16_t,
                                               std::uint_least32_t,
                                               std::uint_least32_t>;

//!\brief Metafunction overload for types that are in seqan3::detail::char_adaptations.
template <typename type_in_list>
//!\cond
    requires meta::in<char_adaptations, type_in_list>::value
//!\endcond
struct is_char_adaptation<type_in_list> :
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

/*!\brief Specialisation of seqan3::underlying_char for char types.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \relates seqan3::char_adaptation_concept
 * \sa seqan3::underlying_char_t
 */
template <typename char_type>
//!\cond
    requires detail::is_char_adaptation_v<char_type>
//!\endcond
struct underlying_char<char_type>
{
    //!\brief The same type as char_type.
    using type = char_type;
};

/*!\brief Specialisation of seqan3::underlying_rank for char types.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \relates seqan3::char_adaptation_concept
 * \sa seqan3::underlying_rank_t
 */
template <typename char_type>
//!\cond
    requires detail::is_char_adaptation_v<char_type>
//!\endcond
struct underlying_rank<char_type>
{
    //!\brief An unsigned integer type of the same size as `char_type`.
    using type = meta::at<detail::char_adaptations_rank_types,
                          meta::find_index<detail::char_adaptations, char_type>>;
};

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

/*!\brief Specialisation of seqan3::alphabet_size that delegates for char types.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \relates seqan3::char_adaptation_concept
 * \sa seqan3::alphabet_size_v
 */
template <typename char_type>
//!\cond
    requires detail::is_char_adaptation_v<char_type>
//!\endcond
struct alphabet_size<char_type>
{
    //!\brief Smallest unsigned integral type that can hold value;
    using type = detail::min_viable_uint_t<static_cast<uint64_t>(std::numeric_limits<char_type>::max()) + 1 -
                                           std::numeric_limits<char_type>::lowest()>;
    //!\brief The alphabet's size.
    static constexpr type value =
        static_cast<type>(std::numeric_limits<char_type>::max()) + 1 - std::numeric_limits<char_type>::lowest();
};

// ------------------------------------------------------------------
// free functions
// ------------------------------------------------------------------

/*!\name Free function wrappers for the char alphabet adaptation
 * \brief For `char`, `char16_t` and `char32_t` do conversion to/from uint types.
 * \ingroup adaptation
 * \relates seqan3::char_adaptation_concept
 * \{
 */
/*!\brief Converting char to char is no-op (it will just return the value you pass in).
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `char`).
 */
template <typename char_type>
constexpr underlying_char_t<char_type> to_char(char_type const chr)
    requires detail::is_char_adaptation_v<char_type>
{
    return chr;
}

/*!\brief Convert char to rank by casting to an unsigned integral type of same size.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename char_type>
constexpr underlying_rank_t<char_type> to_rank(char_type const chr)
    requires detail::is_char_adaptation_v<char_type>
{
    return chr;
}

/*!\brief Assign a char to the char type (same as calling `=`).
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr The alphabet letter that you wish to assign to.
 * \param chr2 The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename char_type>
constexpr char_type & assign_char(char_type & chr, underlying_char_t<char_type> const chr2)
    requires detail::is_char_adaptation_v<char_type>
{
    return chr = chr2;
}

/*!\brief Assign a char to the char type (same as calling `=`).
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr An alphabet letter temporary.
 * \param chr2 The `char` value you wish to assign.
 * \returns The assignment result as a temporary.
 */
template <typename char_type>
constexpr char_type && assign_char(char_type && chr, underlying_char_t<char_type> const chr2)
    requires detail::is_char_adaptation_v<std::remove_reference_t<char_type>>
{
    return std::move(chr = chr2);
}

/*!\brief Assigning a rank to a char is the same as assigning it a numeric value.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr The alphabet letter that you wish to assign to.
 * \param rank The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename char_type>
constexpr char_type & assign_rank(char_type & chr, underlying_rank_t<char_type> const rank)
    requires detail::is_char_adaptation_v<char_type>
{
    return chr = rank;
}

/*!\brief Assigning a rank to a char is the same as assigning it a numeric value.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr An alphabet letter temporary.
 * \param rank The `rank` value you wish to assign.
 * \returns The assignment result as a temporary.
 */
template <typename char_type>
constexpr char_type && assign_rank(char_type && chr, underlying_rank_t<char_type> const rank)
    requires detail::is_char_adaptation_v<std::remove_reference_t<char_type>>
{
    return std::move(chr = rank);
}
//!\}
//!\}

} // namespace seqan3

// concept tests in alphabet/adaptation.hpp because of include order constraints
