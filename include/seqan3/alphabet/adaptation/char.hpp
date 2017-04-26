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

/*!\file alphabet/adaptation/char.hpp
 * \ingroup adaptation
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides alphabet adaptations for standard char types.
 * \details
 * This file provides function and metafunction overloads so that the following types
 * fulfil the seqan3::alphabet_concept:
 *   * `char`
 *   * `wchar_t`
 *   * `char16_t`
 *   * `char32_t`
 *
 * You will likely not use these interfaces directly, they are, however, very helpful
 * for conversions between other alphabets and between other alphabets and characters.
 *
 * \attention
 *   * Note that `signed char` and `unsigned char` are absent from the list, because of
 * their type ambiguity with `int8_t` and `uint8_t`.
 *   * Please be aware that if you also include alphabet/concept.hpp, you need to do so **after**
 * including this file, not before.
 */

#pragma once

#include <limits>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept_fwd.hpp>
#include <seqan3/core/detail/int_types.hpp>

//!\cond
// NOTE(h-2): 'tis some nice code witchery that gives an error if the include order is wrong.
namespace
{

template <std::nullptr_t t>
concept bool alphabet_concept = true;

}

namespace seqan3
{

static_assert(alphabet_concept<nullptr>, "You must include alphabet/concept.hpp AFTER including alphabet/adaptation/char.hpp, not before!");

}
//!\endcond

namespace seqan3::detail
{
//!\addtogroup adaptation
//!\{

//!\brief The list of types that are defined as char adaptations.
using char_adaptations            = meta::list<char,
                                               wchar_t,
                                               char16_t,
                                               char32_t>;
//!\brief The corresponding list of rank types
using char_adaptations_rank_types = meta::list<std::uint_least8_t,
                                               std::uint_least32_t,
                                               std::uint_least16_t,
                                               std::uint_least32_t>;

//!\brief Metafunction that indicates whether a type is a char alphabet adaptation.
template <typename type>
struct is_char_adaptation :
    public std::conditional_t<meta::in<char_adaptations, type>::value, std::true_type, std::false_type>
{};

//!\brief Shortcut for seqan3::detail::is_char_adaptation.
template <typename type>
constexpr bool is_char_adaptation_v = is_char_adaptation<type>::value;

//!\}
} // namespace seqan3::detail

namespace seqan3
{

//!\addtogroup adaptation
//!\{

// ------------------------------------------------------------------
// type metafunctions
// ------------------------------------------------------------------

//!\brief Type metafunction that shall return the type of an alphabet in char representation. [char specialization]
//!\tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
//!\sa seqan3::underlying_char
template <typename char_type>
    requires detail::is_char_adaptation_v<char_type>
struct underlying_char<char_type>
{
    //!\brief The same type as char_type.
    using type = char_type;
};

//!\brief Type metafunction that shall return the type of an alphabet in rank representation. [char specialization]
//!\tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
//!\sa seqan3::underlying_rank
template <typename char_type>
    requires detail::is_char_adaptation_v<char_type>
struct underlying_rank<char_type>
{
    //!\brief An unsigned integer type of the same size as `char_type`.
    using type = meta::at<detail::char_adaptations_rank_types,
                          meta::find_index<detail::char_adaptations, char_type>>;
};

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

//!\brief Value metafunction that shall return the size of an alphabet. [char specialization]
//!\tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
//!\sa seqan3::alphabet_size
template <typename char_type>
    requires detail::is_char_adaptation_v<char_type>
struct alphabet_size<char_type>
{
    //!\brief Smallest unsigned integral type that can hold value;
    using type = detail::min_viable_uint_t<(std::size_t)std::numeric_limits<char_type>::max() + 1 -
                                            std::numeric_limits<char_type>::lowest()>;
    //!\brief The alphabet's size.
    static constexpr type value =
        (type)std::numeric_limits<char_type>::max() + 1 - std::numeric_limits<char_type>::lowest();
};

// ------------------------------------------------------------------
// free functions
// ------------------------------------------------------------------

/*!\name Free function wrappers for the char alphabet adaptation
 * \brief For `char`, `wchar_t`, `char16_t` and `char32_t` do conversion to/from uint types.
 * \ingroup adaptation
 * \{
 */
/*!\brief Converting char to char is no-op (it will just return the value you pass in).
 * \tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
 * \param c The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `char`).
 */
template <typename char_type>
constexpr underlying_char_t<char_type> to_char(char_type const c)
    requires detail::is_char_adaptation_v<char_type>
{
    return c;
}

/*!\brief Convert char to rank by casting to an unsigned integral type of same size.
 * \tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
 * \param c The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename char_type>
constexpr underlying_rank_t<char_type> to_rank(char_type const c)
    requires detail::is_char_adaptation_v<char_type>
{
    return c;
}

/*!\brief The same as calling `=` on the char.
 * \tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
 * \param c The alphabet letter that you wish to assign to.
 * \param in The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename char_type>
constexpr char_type & assign_char(char_type & c, underlying_char_t<char_type> const in)
    requires detail::is_char_adaptation_v<char_type>
{
    return c = in;
}

/*!\brief The same as calling `=` on the char.
 * \tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
 * \param c An alphabet letter temporary.
 * \param in The `char` value you wish to assign.
 * \returns The assignment result as a temporary.
 */
template <typename char_type>
constexpr char_type && assign_char(char_type && c, underlying_char_t<char_type> const in)
    requires detail::is_char_adaptation_v<std::remove_reference_t<char_type>>
{
    return std::move(c = in);
}

/*!\brief Assigning a rank to a char is the same assigning it a numeric value.
 * \tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
 * \param c The alphabet letter that you wish to assign to.
 * \param in The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename char_type>
constexpr char_type & assign_rank(char_type & c, underlying_rank_t<char_type> const in)
    requires detail::is_char_adaptation_v<char_type>
{
    return c = in;
}

/*!\brief Assigning a rank to a char is the same assigning it a numeric value.
 * \tparam char_type One of `char`, `wchar_t`, `char16_t` or `char32_t`.
 * \param c An alphabet letter temporary.
 * \param in The `rank` value you wish to assign.
 * \returns The assignment result as a temporary.
 */
template <typename char_type>
constexpr char_type && assign_rank(char_type && c, underlying_rank_t<char_type> const in)
    requires detail::is_char_adaptation_v<std::remove_reference_t<char_type>>
{
    return std::move(c = in);
}
//!\}
//!\}

} // namespace seqan3

// concept tests in alphabet/adaptation.hpp because of include order constraints
