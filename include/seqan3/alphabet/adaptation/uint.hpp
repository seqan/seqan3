// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
    requires meta::in<uint_adaptations, std::remove_reference_t<type_in_list>>::value
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
                          meta::find_index<detail::uint_adaptations, std::remove_reference_t<uint_type>>>;
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
struct alphabet_size<uint_type> :
    std::integral_constant<detail::min_viable_uint_t<detail::size_in_values_v<uint_type>>,
                           detail::size_in_values_v<uint_type>>
{};

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
constexpr underlying_char_t<uint_type> to_char(uint_type const intgr) noexcept
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
constexpr underlying_rank_t<uint_type> to_rank(uint_type const intgr) noexcept
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
constexpr uint_type & assign_char(uint_type & intgr, underlying_char_t<uint_type> const chr) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr = chr;
}

//!\overload
template <typename uint_type>
constexpr uint_type assign_char(uint_type &&, underlying_char_t<uint_type> const chr) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return chr;
}

//!\brief For adaptations seqan3::assign_char_strict behaves exactly as seqan3::assign_char.
template <typename uint_type>
constexpr uint_type & assign_char_strict(uint_type & intgr, underlying_char_t<uint_type> const chr) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr = chr;
}

//!\overload
template <typename uint_type>
constexpr uint_type assign_char_strict(uint_type &&, underlying_char_t<uint_type> const chr) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return chr;
}

//!\brief For char adaptations, all character values are valid.
template <typename uint_type>
constexpr bool char_is_valid_for(underlying_char_t<uint_type> const) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return true;
}

/*!\brief Assign a rank to to the uint (same as calling `=`).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to assign to.
 * \param intgr2 The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename uint_type>
constexpr uint_type & assign_rank(uint_type & intgr, underlying_rank_t<uint_type> const intgr2) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr = intgr2;
}

//!\overload
template <typename uint_type>
constexpr uint_type assign_rank(uint_type &&, underlying_rank_t<uint_type> const intgr2) noexcept
    requires detail::is_uint_adaptation_v<uint_type>
{
    return intgr2;
}
//!\}
//!\}

} // namespace seqan3

// concept tests in alphabet/adaptation.hpp because of include order constraints
