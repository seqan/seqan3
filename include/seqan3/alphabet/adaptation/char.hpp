// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
    requires meta::in<char_adaptations, std::remove_reference_t<type_in_list>>::value
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
                          meta::find_index<detail::char_adaptations, std::remove_reference_t<char_type>>>;
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
struct alphabet_size<char_type> :
    std::integral_constant<detail::min_viable_uint_t<detail::size_in_values_v<char_type>>,
                           detail::size_in_values_v<char_type>>
{};

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
constexpr underlying_char_t<char_type> to_char(char_type const chr) noexcept
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
constexpr underlying_rank_t<char_type> to_rank(char_type const chr) noexcept
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
constexpr char_type & assign_char(char_type & chr, underlying_char_t<char_type> const chr2) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return chr = chr2;
}

//!\overload
template <typename char_type>
constexpr char_type assign_char(char_type &&, underlying_char_t<char_type> const chr2) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return chr2;
}

/*!\brief For adaptations seqan3::assign_char_strict behaves exactly as seqan3::assign_char.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr The alphabet letter that you wish to assign to.
 * \param chr2 The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename char_type>
constexpr char_type & assign_char_strict(char_type & chr, underlying_char_t<char_type> const chr2) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return assign_char(chr, chr2);
}

//!\overload
template <typename char_type>
constexpr char_type assign_char_strict(char_type &&, underlying_char_t<char_type> const chr2) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return assign_char(char_type{}, chr2);
}

//!\brief For char adaptations, all character values are valid.
template <typename char_type>
constexpr bool char_is_valid_for(underlying_char_t<char_type> const) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return true;
}

/*!\brief Assigning a rank to a char is the same as assigning it a numeric value.
 * \tparam char_type One of `char`, `char16_t`, `char32_t` or `wchar_t`.
 * \param chr The alphabet letter that you wish to assign to.
 * \param rank The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename char_type>
constexpr char_type & assign_rank(char_type & chr, underlying_rank_t<char_type> const rank) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return chr = rank;
}

//!\overload
template <typename char_type>
constexpr char_type assign_rank(char_type &&, underlying_rank_t<char_type> const rank) noexcept
    requires detail::is_char_adaptation_v<char_type>
{
    return rank;
}
//!\}
//!\}

} // namespace seqan3

// concept tests in alphabet/adaptation.hpp because of include order constraints
