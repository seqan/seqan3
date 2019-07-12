// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides alphabet adaptations for standard uint types.
 * \details
 * This file provides function and type trait overloads so that the following types
 * fulfil the seqan3::Alphabet:
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

#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{
//!\addtogroup adaptation
//!\{

//!\brief Whether a type is `uint8_t`, `uint16_t` or `uint32_t`.
template <typename type, typename type2 = std::remove_reference_t<type>>
constexpr bool is_uint_adaptation_v = std::Same<type2, uint8_t>  ||
                                      std::Same<type2, uint16_t> ||
                                      std::Same<type2, uint32_t>;
//!\}
} // namespace seqan3::detail

namespace seqan3::custom
{

/*!\name Free function wrappers for the uint alphabet adaptation
 * \brief For `uint8_t`, `uint16_t` and `uint32_t` do conversion to/from char types.
 * \ingroup adaptation
 * \{
 */

/*!\brief Return the number of values the uint type can take.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The parameter's actual value is ignored.
 * \returns The respective size (e.g. 256 for `uint8_t`).
 */
template <typename uint_type>
//!\cond
    requires ::seqan3::detail::is_uint_adaptation_v<uint_type>
//!\endcond
constexpr auto alphabet_size(uint_type const & SEQAN3_DOXYGEN_ONLY(intgr)) noexcept
{
    return detail::min_viable_uint_t<detail::size_in_values_v<uint_type>>{detail::size_in_values_v<uint_type>};
}

/*!\brief Converting uint to char casts to a character type of same size.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `uint`).
 */
template <typename uint_type>
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
constexpr auto to_char(uint_type const intgr) noexcept
{
    if constexpr (std::Same<uint_type, uint8_t>)
        return static_cast<char>(intgr);
    else if constexpr (std::Same<uint_type, uint16_t>)
        return static_cast<char16_t>(intgr);
    else
        return static_cast<char32_t>(intgr);
}

/*!\brief Converting uint to rank is a no-op (it will just return the value you pass in).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to convert to rank.
 * \returns `intgr`.
 */
template <typename uint_type>
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
constexpr uint_type to_rank(uint_type const intgr) noexcept
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
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
constexpr uint_type & assign_char_to(decltype(to_char(uint_type{})) const chr, uint_type & intgr) noexcept
{
    return intgr = chr;
}

/*!\brief Assign a rank to to the uint (same as calling `=`).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param intgr The alphabet letter that you wish to assign to.
 * \param intgr2 The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename uint_type>
//!\cond
    requires detail::is_uint_adaptation_v<uint_type>
//!\endcond
constexpr uint_type & assign_rank_to(uint_type const intgr2, uint_type & intgr) noexcept
{
    return intgr = intgr2;
}

//!\}
} // namespace seqan3::custom
